// implementation for ds_basis

#include <cassert>
#include <cmath>

#include <numerics/matrix_decomp.h>

#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/dm_mask.h>
#include <Rd/refinable.h>
#include <Rd/r_index.h>

namespace WaveletTL
{
  template <int d, int dT>
  DSBasis<d,dT>::DSBasis(const int s0, const int s1, const int sT0, const int sT1,
			 DSBiorthogonalizationMethod bio) {
    assert(std::max(s0,s1) < d && std::max(sT0,sT1) < dT);
    
    assert(((s0 == 0 && sT0 > 0) || (s0 > 0 && sT0 == 0))
	   && ((s1 == 0 && sT1 > 0) || (s1 > 0 && sT1 == 0)));
    
    this->s0 = s0;
    this->s1 = s1;
    this->sT0 = sT0;
    this->sT1 = sT1;
    this->bio = bio;

    setup_GammaLR();
    setup_CX_CXT();
    setup_CXA_CXAT();

    // IGPMlib reference: I_Basis_Bspline_s::Setup()

    setup_Cj();

    Matrix<double> ml = ML(); // (3.5.2)
    Matrix<double> mr = MR(); // (3.5.2)
    SparseMatrix<double> mj0;   setup_Mj0  (ml,   mr,   mj0);   // (3.5.1)

    Matrix<double> mltp = MLTp(); // (3.5.6)
    Matrix<double> mrtp = MRTp(); // (3.5.6)
    SparseMatrix<double> mj0tp; setup_Mj0Tp(mltp, mrtp, mj0tp); // (3.5.5)
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::alpha(const int m, const unsigned int r) const {
    double result = 0;
    if (r == 0)
      return 1; // [DKU] (5.1.1)
    else {
      if (m == 0) {
	// [DKU] (5.1.3)
	for (int k = ell1<d>(); k <= ell2<d>(); k++) {
	  double help = 0;
	  for (unsigned int s = 0; s < r; s++)
	    help += binomial(r, s) * intpower(k, r-s) * alpha(0, s);
	  result += cdf.a().get_coefficient(MultiIndex<int,1>(k)) * help;
	}
	result /= (double)((1<<(r+1))-2);
      } else {
	// [DKU] (5.1.2)
	for (unsigned int i = 0; i <= r; i++)
	  result += binomial(r, i) * intpower(m, i) * alpha(0, r-i);
      }
    }
    return result;
  }
  
  template <int d, int dT>
  const double
  DSBasis<d,dT>::alphaT(const int m, const unsigned int r) const {
    double result = 0;
    if (r == 0)
      return 1; // [DKU] (5.1.1)
    else {
      if (m == 0) {
	// [DKU] (5.1.3)
	for (int k = ell1T<d,dT>(); k <= ell2T<d,dT>(); k++) {
	  double help = 0;
	  for (unsigned int s = 0; s < r; s++)
	    help += binomial(r, s) * intpower(k, r-s) * alphaT(0, s);
	  result += cdf.aT().get_coefficient(MultiIndex<int,1>(k)) * help;
	}
	result /= (double)((1<<(r+1))-2);
      } else {
	// [DKU] (5.1.2)
	for (unsigned int i = 0; i <= r; i++)
	  result += binomial(r, i) * intpower(m, i) * alphaT(0, r-i);
      }
    }
    return result;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaL(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_l(); q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaLT(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2<d>())/2.0); q < ell_l(); q++)
      result += alphaT(q, r) * cdf.a().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaR(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2T<d,dT>())/2.0); q < ellT_r(); q++)
      result += alpha(q, r) * cdf.aT().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  const double
  DSBasis<d,dT>::betaRT(const int m, const unsigned int r) const {
    // [DKU] (3.2.31)
    double result = 0;
    for (int q = (int)ceil((m-ell2<d>())/2.0); q < ell_r(); q++)
      result += alphaT(q, r) * cdf.a().get_coefficient(MultiIndex<int,1>(m-2*q));
    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_GammaLR() {
    // IGPMlib reference: I_Mask_Bspline::EvalGammaL(), ::EvalGammaR()

    const unsigned int GammaLsize = std::max(d-s0, dT-sT0);
    GammaL.resize(GammaLsize, GammaLsize);

    const unsigned int GammaRsize = std::max(d-s1, dT-sT1);
    GammaR.resize(GammaRsize, GammaRsize);

    // 1. compute the integrals
    //      z(s,t) = \int_0^1\phi(x-s)\tilde\phi(x-t)\,dx
    //    exactly with the [DM] trick
    MultivariateRefinableFunction<DMMask2<HaarMask, CDFMask_primal<d>, CDFMask_dual<d,dT> >,2> zmask;
    InfiniteVector<double, MultiIndex<int,2> > zvalues(zmask.evaluate());

    // 2. compute the integrals
    //      I(nu,mu) = \int_0^\infty\phi(x-\nu)\tilde\phi(x-\mu)\,dx
    //    exactly using the z(s,t) values

    const int I1Low  = -ell2<d>()+1;
    const int I1Up   = ell_l()-d+dT-1;
    const int I1Up_r = ell_r()-d+dT-1;
    const int I2Low  = -ell2T<d,dT>()+1;
    const int I2Up   = ellT_l()-1;
    const int I2Up_r = ellT_r()-1;

    Matrix<double> I(std::max(I1Up,I1Up_r)-I1Low+1,
		     std::max(I2Up,I2Up_r)-I2Low+1);
    for (int nu = I1Low; nu < -ell1<d>(); nu++)
      for (int mu = I2Low; mu < -ell1T<d,dT>(); mu++) {
  	double help(0);
  	const int diff(mu - nu);
   	const int sLow = std::max(-ell2<d>()+1,-ell2T<d,dT>()+1-diff);
  	for (int s = sLow; s <= nu; s++)
  	  help += zvalues.get_coefficient(MultiIndex<int,2>(s, s+diff));
  	I(-I1Low+nu, -I2Low+mu) = help; // [DKU] (5.1.7)
      }
    for (int nu = I1Low; nu <= std::max(I1Up,I1Up_r); nu++)
      for (int mu = I2Low; mu <= std::max(I2Up,I2Up_r); mu++) {
 	if ((nu >= -ell1<d>()) || ((nu <= ell_l()-1) && (mu >= -ell1T<d,dT>())))
 	  I(-I1Low+nu, -I2Low+mu) = (nu == mu ? 1 : 0); // [DKU] (5.1.6)
      }

//     cout << "I=" << endl << I << endl;

    // 3. finally, compute the Gramian GammaL
    for (int r = s0; r < d; r++) {
      for (int k = sT0; k < dT; k++) {
     	double help = 0;
     	for (int nu = I1Low; nu < ell_l(); nu++)
     	  for (int mu = I2Low; mu <= I2Up; mu++)
     	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
    	GammaL(r-s0, k-sT0) = help; // [DKU] (5.1.4)
      }
      for (int k = dT; k < d+sT0-s0; k++) {
	double help = 0;
	for (int nu = I1Low; nu < ell_l(); nu++)
	  help += alphaT(nu, r) * I(-I1Low+nu, -I2Low+ellT_l()-dT-s0-sT0+k);
	GammaL(r-s0, k-sT0) = help; // analogous to [DKU] (5.1.4) (!)
      }
    }
    
    for (int r = d; r < dT+s0-sT0; r++)
      for (int k = sT0; k < dT; k++) {
    	double help = 0;
    	for (int mu = I2Low; mu <= I2Up; mu++)
    	  help += alphaT(mu, k) * I(-I1Low+ell_l()-d-s0-sT0+r, -I2Low+mu);
    	GammaL(r-s0, k-sT0) = help; // [DKU] (5.1.5)
      }

//     cout << "GammaL=" << endl << GammaL << endl;

#if 0
    Matrix<double> Gammafull(dT);
    for (int r = 0; r < d; r++)
      for (int k = 0; k < dT; k++) {
 	double help = 0;
 	for (int nu = I1Low; nu < ell_l(); nu++)
 	  for (int mu = I2Low; mu <= I2Up; mu++)
 	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
 	Gammafull(r, k) = help;
      }
    cout << "For testing, the full Gramian (without b.c.'s):" << endl
 	 << Gammafull;
#endif

    // The same for GammaR:

    for (int r = s1; r < d; r++) {
      for (int k = sT1; k < dT; k++) {
     	double help = 0;
     	for (int nu = I1Low; nu < ell_r(); nu++)
     	  for (int mu = I2Low; mu <= I2Up_r; mu++)
     	    help += alphaT(nu, r) * alpha(mu, k) * I(-I1Low+nu, -I2Low+mu);
    	GammaR(r-s1, k-sT1) = help; // [DKU] (5.1.4)
      }
      for (int k = dT; k < d+sT1-s1; k++) {
	double help = 0;
	for (int nu = I1Low; nu < ell_r(); nu++)
	  help += alphaT(nu, r) * I(-I1Low+nu, -I2Low+ellT_r()-dT-s1-sT1+k);
	GammaR(r-s1, k-sT1) = help; // analogous to [DKU] (5.1.4) (!)
      }
    }
    
    for (int r = d; r < dT+s1-sT1; r++)
      for (int k = sT1; k < dT; k++) {
    	double help = 0;
    	for (int mu = I2Low; mu <= I2Up_r; mu++)
    	  help += alphaT(mu, k) * I(-I1Low+ell_r()-d-s1-sT1+r, -I2Low+mu);
    	GammaR(r-s1, k-sT1) = help; // [DKU] (5.1.5)
      }

//     cout << "GammaR=" << endl << GammaR << endl;
  }

  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_CX_CXT()
  {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    CL.resize(GammaL.row_dimension(), GammaL.column_dimension());
    CLT.resize(GammaL.row_dimension(), GammaL.column_dimension());
    CR.resize(GammaR.row_dimension(), GammaR.column_dimension());
    CRT.resize(GammaR.row_dimension(), GammaR.column_dimension());

    if (bio == none) {
      CL.diagonal(GammaL.row_dimension(), 1.0);
      Matrix<double> CLGammaLInv;
      QUDecomposition<double>(GammaL).inverse(CLGammaLInv);
      CLT = transpose(CLGammaLInv);

      CR.diagonal(GammaR.row_dimension(), 1.0);
      Matrix<double> CRGammaRInv;
      QUDecomposition<double>(GammaR).inverse(CRGammaRInv);
      CRT = transpose(CRGammaRInv);      
    }
    
    if (bio == SVD) {
      MathTL::SVD<double> svd(GammaL);
      Matrix<double> U, V;
      Vector<double> S;
      svd.getU(U);
      svd.getV(V);
      svd.getS(S);
      for (unsigned int i = 0; i < GammaL.row_dimension(); i++) {
 	S[i] = 1.0 / sqrt(S[i]);
 	for (unsigned int j = 0; j < GammaL.row_dimension(); j++)
 	  CL(i, j)  = S[i] * U(j, i);
      }
      CL.compress(1e-14);
      Matrix<double> CLGammaLInv;
      QUDecomposition<double>(CL * GammaL).inverse(CLGammaLInv);
      CLT = transpose(CLGammaLInv);
      CLT.compress(1e-14);
      
      MathTL::SVD<double> svd_r(GammaR);
      svd_r.getU(U);
      svd_r.getV(V);
      svd_r.getS(S);
      for (unsigned int i = 0; i < GammaR.row_dimension(); i++) {
 	S[i] = 1.0 / sqrt(S[i]);
 	for (unsigned int j = 0; j < GammaR.row_dimension(); j++)
 	  CR(i, j)  = S[i] * U(j, i);
      }
      CR.compress(1e-14);
      Matrix<double> CRGammaRInv;
      QUDecomposition<double>(CR * GammaR).inverse(CRGammaRInv);
      CRT = transpose(CRGammaRInv);      
      CRT.compress(1e-14);
    }
    
//     if (bio == Bernstein) {
//       double b(0);
//       if (d == 2)
//  	b = 0.7; // cf. [DKU]
//       else {
//  	b = 3.6;
//       }
      
//       // CL : transformation to Bernstein polynomials (triangular)
//       //   (Z_b^m)_{r,l} = (-1)^{r-1}\binom{m-1}{r}\binom{r}{l}b^{-r}, r>=l, 0 otherwise
//       for (int j(0); j < d; j++)
//  	for (int k(0); k <= j; k++)
//  	  CL_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < lupT-llowT+1; j++)
//  	CL_(j, j) = 1.0;
      
//       Matrix<double> CLGammaLInv;
//       QUDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
//       CLT_ = transpose(CLGammaLInv);
      
//       CLT_.compress(1e-14);

//       for (int j(0); j < d; j++)
// 	for (int k(0); k <= j; k++)
// 	  CR_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < rlowTh-rupTh+1; j++)
// 	CR_(j, j) = 1.0;
      
//       Matrix<double> CRGammaRInv;
//       QUDecomposition<double>(CR_ * GammaR_).inverse(CRGammaRInv);
//       CRT_ = transpose(CRGammaRInv);
      
//       CRT_.compress(1e-14);
//     }
    
    if (bio == partialSVD) {
      MathTL::SVD<double> svd(GammaL);
      Matrix<double> U, V;
      Vector<double> S;
      svd.getU(U);
      svd.getV(V);
      svd.getS(S);
      Matrix<double> GammaLInv;
      QUDecomposition<double>(GammaL).inverse(GammaLInv);
      const double a = 1.0 / GammaLInv(0, 0);
      Matrix<double> R(GammaL.row_dimension(), GammaL.column_dimension());
      for (unsigned int i = 0; i < R.row_dimension(); i++)
 	S[i] = sqrt(S[i]);
      for (unsigned int j = 0; j < R.row_dimension(); j++)
 	R(0, j) = a * V(0, j) / S[j];
      for (unsigned int i = 1; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  R(i, j) = U(j, i) * S[j];
      for (unsigned int i = 0; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  U(i, j) /= S[i];
      CL = R*U;
      CL.compress(1e-14);
      Matrix<double> CLGammaLInv;
      QUDecomposition<double>(CL * GammaL).inverse(CLGammaLInv);
      CLT = transpose(CLGammaLInv);
      CLT.compress(1e-14);

      MathTL::SVD<double> svd_r(GammaR);
      svd_r.getU(U);
      svd_r.getV(V);
      svd_r.getS(S);
      Matrix<double> GammaRInv;
      QUDecomposition<double>(GammaR).inverse(GammaRInv);
      const double a_r = 1.0 / GammaRInv(0, 0);
      R.resize(GammaR.row_dimension(), GammaR.column_dimension());
      for (unsigned int i = 0; i < R.row_dimension(); i++)
	S[i] = sqrt(S[i]);
      for (unsigned int j = 0; j < R.row_dimension(); j++)
 	R(0, j) = a_r * V(0, j) / S[j];
      for (unsigned int i = 1; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  R(i, j) = U(j, i) * S[j];
      for (unsigned int i = 0; i < R.row_dimension(); i++)
 	for (unsigned int j = 0; j < R.row_dimension(); j++)
 	  U(i, j) /= S[i];
      CR = R*U;
      CR.compress(1e-14);
      Matrix<double> CRGammaRInv;
      QUDecomposition<double>(CR * GammaR).inverse(CRGammaRInv);
      CRT = transpose(CRGammaRInv);
      CRT.compress(1e-14);
    }

//     if (bio_ == BernsteinSVD) {
//       double b(0);
//       if (d == 2)
// 	b = 0.7; // cf. [DKU]
//       else {
// 	b = 3.6;
//       }
      
//       // CL_ : transformation to Bernstein polynomials (triangular)
//       //   (Z_b^m)_{r,l} = (-1)^{r-1}\binom{m-1}{r}\binom{r}{l}b^{-r}, r>=l, 0 otherwise
//       for (int j(0); j < d; j++)
// 	for (int k(0); k <= j; k++)
// 	  CL_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < lupT-llowT+1; j++)
// 	CL_(j, j) = 1.0;
      
//       Matrix<double> GammaLNew(CL_ * GammaL_);
      
//       MathTL::SVD<double> svd(GammaLNew);
//       Matrix<double> U, V;
//       Vector<double> S;
//       svd.getU(U);
//       svd.getV(V);
//       svd.getS(S);
      
//       Matrix<double> GammaLNewInv;
//       QUDecomposition<double>(GammaLNew).inverse(GammaLNewInv);
//       const double a = 1.0 / GammaLNewInv(0, 0);
      
//       Matrix<double> R(lupT-llowT+1, lupT-llowT+1);
//       for (int i(0); i < lupT-llowT+1; i++)
// 	S[i] = sqrt(S[i]);
//       for (int j(0); j < lupT-llowT+1; j++)
// 	R(0, j) = a * V(0, j) / S[j];
//       for (int i(1); i < lupT-llowT+1; i++)
// 	for (int j(0); j < lupT-llowT+1; j++)
// 	  R(i, j) = U(j, i) * S[j];
      
//       for (int i(0); i < lupT-llowT+1; i++)
// 	for (int j(0); j < lupT-llowT+1; j++)
// 	  U(i, j) /= S[i];
//       CL_ = R*U*CL_;
      
//       CL_.compress(1e-14);

//       Matrix<double> CLGammaLInv;
//       QUDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
//       CLT_ = transpose(CLGammaLInv);

//       CLT_.compress(1e-14);

//       for (int j(0); j < d; j++)
// 	for (int k(0); k <= j; k++)
// 	  CR_(k, j) = minus1power(j-k) * std::pow(b, -j) * binomial(d-1, j) * binomial(j, k);
//       for (int j(d); j < rlowTh-rupTh+1; j++)
// 	CR_(j, j) = 1.0;

//       Matrix<double> GammaRNew(CR_ * GammaR_);
      
//       MathTL::SVD<double> svd_r(GammaRNew);
//       svd_r.getU(U);
//       svd_r.getV(V);
//       svd_r.getS(S);

//       Matrix<double> GammaRNewInv;
//       QUDecomposition<double>(GammaRNew).inverse(GammaRNewInv);
//       const double a_r = 1.0 / GammaRNewInv(0, 0);
      
//       R.resize(rlowTh-rupTh+1,rlowTh-rupTh+1);
//       for (int i(0); i < rlowTh-rupTh+1; i++)
// 	S[i] = sqrt(S[i]);
//       for (int j(0); j < rlowTh-rupTh+1; j++)
// 	R(0, j) = a_r * V(0, j) / S[j];
//       for (int i(1); i < rlowTh-rupTh+1; i++)
// 	for (int j(0); j < rlowTh-rupTh+1; j++)
// 	  R(i, j) = U(j, i) * S[j];
      
//       for (int i(0); i < rlowTh-rupTh+1; i++)
// 	for (int j(0); j < rlowTh-rupTh+1; j++)
// 	  U(i, j) /= S[i];
//       CR_ = R*U*CR_;

//       CR_.compress(1e-14);

//       Matrix<double> CRGammaRInv;
//       QUDecomposition<double>(CR_ * GammaR_).inverse(CRGammaRInv);
//       CRT_ = transpose(CRGammaRInv);

//       CRT_.compress(1e-14);
//     }
    
#if 0
    // check biorthogonality of the matrix product CL * GammaL * (CLT)^T
    cout << "GammaL=" << endl << GammaL << endl;
    cout << "CL=" << endl << CL << endl;
    cout << "CLT=" << endl << CLT << endl;
    Matrix<double> check(CL * GammaL * transpose(CLT));
    for (unsigned int i(0); i < check.row_dimension(); i++)
      check(i, i) -= 1;
    cout << "error for CLT: " << row_sum_norm(check) << endl;
#endif

    QUDecomposition<double>(CL).inverse(inv_CL);
    QUDecomposition<double>(CLT).inverse(inv_CLT);

#if 0
    // check biorthogonality of the matrix product CR * GammaR * (CRT)^T
    cout << "GammaR=" << endl << GammaR << endl;
    cout << "CR=" << endl << CR << endl;
    cout << "CRT=" << endl << CRT << endl;
    Matrix<double> check2(CR * GammaR * transpose(CRT));
    for (unsigned int i(0); i < check2.row_dimension(); i++)
      check2(i, i) -= 1;
    cout << "error for CRT: " << row_sum_norm(check2) << endl;
#endif

    QUDecomposition<double>(CR).inverse(inv_CR);
    QUDecomposition<double>(CRT).inverse(inv_CRT);
  }

  template <int d, int dT>
  void DSBasis<d,dT>::setup_CXA_CXAT() {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    // setup CLA <-> AlphaT * (CL)^T
    CLA.resize(ell_l()+ell2<d>()+std::max(0,dT-d+s0-sT0)-1, CL.row_dimension());
    for (int i = 1-ell2<d>(); i < ell_l(); i++)
      for (unsigned int r = s0; r < d; r++)
	CLA(ell2<d>()-1+i, r-s0) = alphaT(i, r);
    for (int i = ell_l(); i < ell_l()+std::max(0,dT-d+s0-sT0); i++)
      CLA(ell2<d>()-1+i, i-ell_l()+d-s0) = 1.0;

//     cout << "CLA before biorthogonalization:" << endl << CLA << endl;
    CLA = CLA * transpose(CL);
    CLA.compress(1e-12);

    // setup CLAT <-> Alpha * (CLT)^T
    CLAT.resize(ellT_l()+ell2T<d,dT>()+std::max(0,d-dT+sT0-s0)-1, CLT.row_dimension());
    for (int i = 1-ell2T<d,dT>(); i < ellT_l(); i++)
      for (unsigned int r = sT0; r < dT; r++)
	CLAT(ell2T<d,dT>()-1+i, r-sT0) = alpha(i, r);
    for (int i = ellT_l(); i < ellT_l()+std::max(0,d-dT+sT0-s0); i++)
      CLAT(ell2T<d,dT>()-1+i, i-ellT_l()+dT-sT0) = 1.0;

//     cout << "CLAT before biorthogonalization:" << endl << CLAT << endl;
    CLAT = CLAT * transpose(CLT);
    CLAT.compress(1e-12);

    // the same for CRA, CRAT:
    CRA.resize(ell_r()+ell2<d>()+std::max(0,dT-d+s1-sT1)-1, CR.row_dimension());
    for (int i = 1-ell2<d>(); i < ell_r(); i++)
      for (unsigned int r = s1; r < d; r++)
	CRA(ell2<d>()-1+i, r-s1) = alphaT(i, r);
    for (int i = ell_r(); i < ell_r()+std::max(0,dT-d+s1-sT1); i++)
      CRA(ell2<d>()-1+i, i-ell_r()+d-s1) = 1.0;

//     cout << "CRA before biorthogonalization:" << endl << CRA << endl;
    CRA = CRA * transpose(CR);
    CRA.compress(1e-12);

    CRAT.resize(ellT_r()+ell2T<d,dT>()+std::max(0,d-dT+sT1-s1)-1, CRT.row_dimension());
    for (int i = 1-ell2T<d,dT>(); i < ellT_r(); i++)
      for (unsigned int r = sT1; r < dT; r++)
	CRAT(ell2T<d,dT>()-1+i, r-sT1) = alpha(i, r);
    for (int i = ellT_r(); i < ellT_r()+std::max(0,d-dT+sT1-s1); i++)
      CRAT(ell2T<d,dT>()-1+i, i-ellT_r()+dT-sT1) = 1.0;

//     cout << "CRAT before biorthogonalization:" << endl << CRAT << endl;
    CRAT = CRAT * transpose(CRT);
    CRAT.compress(1e-12);

#if 0
    cout << "CLA=" << endl << CLA << endl;
    cout << "CLAT=" << endl << CLAT << endl;
    cout << "CRA=" << endl << CRA << endl;
    cout << "CRAT=" << endl << CRAT << endl;
#endif
  }

  template <int d, int dT>
  void DSBasis<d,dT>::setup_Cj() {
    // IGPMlib reference: I_Basis_Bspline_s::setup_Cj(), ::put_Mat()

    // (5.2.5)
    Cj.diagonal(Deltasize(j0()), 1.0);
    Cj.set_block(0, 0, CL);
    Cj.set_block(Deltasize(j0())-CR.row_dimension(),
		 Deltasize(j0())-CR.column_dimension(),
		 CR, true);

    inv_Cj.diagonal(Deltasize(j0()), 1.0);
    inv_Cj.set_block(0, 0, inv_CL);
    inv_Cj.set_block(Deltasize(j0())-inv_CR.row_dimension(),
		     Deltasize(j0())-inv_CR.column_dimension(),
		     inv_CR, true);

    CjT.diagonal(Deltasize(j0()), 1.0);
    CjT.set_block(0, 0, CLT);
    CjT.set_block(Deltasize(j0())-CRT.row_dimension(),
		   Deltasize(j0())-CRT.column_dimension(),
		   CRT, true);

    inv_CjT.diagonal(Deltasize(j0()), 1.0);
    inv_CjT.set_block(0, 0, inv_CLT);
    inv_CjT.set_block(Deltasize(j0())-inv_CRT.row_dimension(),
		      Deltasize(j0())-inv_CRT.column_dimension(),
		      inv_CRT, true);

    Cjp.diagonal(Deltasize(j0()+1), 1.0);
    Cjp.set_block(0, 0, CL);
    Cjp.set_block(Deltasize(j0()+1)-CR.row_dimension(),
		  Deltasize(j0()+1)-CR.column_dimension(),
		  CR, true);

    inv_Cjp.diagonal(Deltasize(j0()+1), 1.0);
    inv_Cjp.set_block(0, 0, inv_CL);
    inv_Cjp.set_block(Deltasize(j0()+1)-inv_CR.row_dimension(),
		      Deltasize(j0()+1)-inv_CR.column_dimension(),
		      inv_CR, true);

    CjpT.diagonal(Deltasize(j0()+1), 1.0);
    CjpT.set_block(0, 0, CLT);
    CjpT.set_block(Deltasize(j0()+1)-CRT.row_dimension(),
		   Deltasize(j0()+1)-CRT.column_dimension(),
		   CRT, true);

    inv_CjpT.diagonal(Deltasize(j0()+1), 1.0);
    inv_CjpT.set_block(0, 0, inv_CLT);
    inv_CjpT.set_block(Deltasize(j0()+1)-inv_CRT.row_dimension(),
		       Deltasize(j0()+1)-inv_CRT.column_dimension(),
		       inv_CRT, true);

#if 0
    cout << "DSBasis: testing setup of Cj:" << endl;

    SparseMatrix<double> test1 = CjT * inv_CjT;
    for (unsigned int i = 0; i < test1.row_dimension(); i++)
      test1.set_entry(i, i, test1.get_entry(i, i) - 1.0);
    cout << "* ||CjT*inv_CjT-I||_1: " << column_sum_norm(test1) << endl;
    cout << "* ||CjT*inv_CjT-I||_infty: " << row_sum_norm(test1) << endl;

    SparseMatrix<double> test2 = Cj * inv_Cj;
    for (unsigned int i = 0; i < test2.row_dimension(); i++)
      test2.set_entry(i, i, test2.get_entry(i, i) - 1.0);
    cout << "* ||Cj*inv_Cj-I||_1: " << column_sum_norm(test2) << endl;
    cout << "* ||Cj*inv_Cj-I||_infty: " << row_sum_norm(test2) << endl;

    SparseMatrix<double> test3 = CjpT * inv_CjpT;
    for (unsigned int i = 0; i < test3.row_dimension(); i++)
      test3.set_entry(i, i, test3.get_entry(i, i) - 1.0);
    cout << "* ||CjpT*inv_CjpT-I||_1: " << column_sum_norm(test3) << endl;
    cout << "* ||CjpT*inv_CjpT-I||_infty: " << row_sum_norm(test3) << endl;

    SparseMatrix<double> test4 = Cjp * inv_Cjp;
    for (unsigned int i = 0; i < test4.row_dimension(); i++)
      test4.set_entry(i, i, test4.get_entry(i, i) - 1.0);
    cout << "* ||Cjp*inv_Cjp-I||_1: " << column_sum_norm(test4) << endl;
    cout << "* ||Cjp*inv_Cjp-I||_infty: " << row_sum_norm(test4) << endl;
#endif
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::ML() const {
    // IGPMlib reference: I_Basis_Bspline_s::ML()
    
    Matrix<double> ML(d+ell_l()+ell2<d>()-1-s0, d-s0);

    for (int row = s0; row < d; row++)
      ML(row-s0, row-s0) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ell_l(); m <= 2*ell_l()+ell1<d>()-1; m++)
      for (int k = s0; k < d; k++)
     	ML(-ell_l()+d+m-s0, k-s0) = alphaT(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*ell_l()+ell1<d>(); m <= ell2<d>()+2*ell_l()-2; m++)
      for (int k = s0; k < d; k++)
	ML(-ell_l()+d+m-s0, k-s0) = betaLT(m, k);

    cout << "ML=" << endl << ML << endl;
    
    return ML;
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::MR() const {
    // IGPMlib reference: I_Basis_Bspline_s::MR()

    Matrix<double> MR(d+ell_r()+ell2<d>()-1-s1, d-s1);

    for (int row = s1; row < d; row++)
      MR(row-s1, row-s1) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ell_r(); m <= 2*ell_r()+ell1<d>()-1; m++)
      for (int k = s1; k < d; k++)
     	MR(-ell_r()+d+m-s1, k-s1) = alphaT(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*ell_r()+ell1<d>(); m <= ell2<d>()+2*ell_r()-2; m++)
      for (int k = s1; k < d; k++)
	MR(-ell_r()+d+m-s1, k-s1) = betaRT(m, k);

    cout << "MR=" << endl << MR << endl;

    return MR;
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::MLTp() const {
    // IGPMlib reference: I_Basis_Bspline_s::MLts()

    Matrix<double> MLTp(dT+ellT_l()+ell2T<d,dT>()-1-sT0, dT-sT0);

    for (int row = sT0; row < dT; row++)
      MLTp(row-sT0, row-sT0) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ellT_l(); m <= 2*ellT_l()+ell1T<d,dT>()-1; m++)
      for (int k = sT0; k < dT; k++)
      	MLTp(-ellT_l()+dT+m-sT0, k-sT0) = alpha(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*ellT_l()+ell1T<d,dT>(); m <= ell2T<d,dT>()+2*ellT_l()-2; m++)
      for (int k = sT0; k < dT; k++)
 	MLTp(-ellT_l()+dT+m-sT0, k-sT0) = betaL(m, k);

    cout << "MLTp=" << endl << MLTp << endl;
    
    return MLTp;
  }

  template <int d, int dT>
  Matrix<double>
  DSBasis<d, dT>::MRTp() const {
    // IGPMlib reference: I_Basis_Bspline_s::MRts()

    Matrix<double> MRTp(dT+ellT_r()+ell2T<d,dT>()-1-sT1, dT-sT1);

    for (int row = sT1; row < dT; row++)
      MRTp(row-sT1, row-sT1) = 1.0 / sqrt(ldexp(1.0, 2*row+1));
    for (int m = ellT_r(); m <= 2*ellT_r()+ell1T<d,dT>()-1; m++)
      for (int k = sT1; k < dT; k++)
     	MRTp(-ellT_r()+dT+m-sT1, k-sT1) = alpha(m, k) / sqrt(ldexp(1.0, 2*k+1));
    for (int m = 2*ellT_r()+ell1T<d,dT>(); m <= ell2T<d,dT>()+2*ellT_r()-2; m++)
      for (int k = sT1; k < dT; k++)
	MRTp(-ellT_r()+dT+m-sT1, k-sT1) = betaR(m, k);
  
    cout << "MRTp=" << endl << MRTp << endl;
 
    return MRTp;
  }

  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_Mj0(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& Mj0) {
    // IGPMlib reference: I_Basis_Bspline_s::Mj0()
    
    // TODO: enhance readability! (<-> [DKU section 3.5])

    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);
    Mj0.resize(njp, nj);

    for (unsigned int i = 0; i < ML.row_dimension(); i++)
      for (unsigned int k = 0; k < ML.column_dimension(); k++)
 	Mj0.set_entry(i, k, ML.get_entry(i, k));
    
    for (unsigned int i = 0; i < MR.row_dimension(); i++)
      for (unsigned int k = 0; k < MR.column_dimension(); k++)
 	Mj0.set_entry(njp-i-1, nj-k-1, MR.get_entry(i, k));
    
    int startrow = d+ell_l()+ell1<d>()-2*s0;
    for (int col = d-s0; col < nj-(d-s1); col++, startrow+=2)
      {
	int row = startrow;
	for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf.a().begin());
	     it != cdf.a().end(); ++it, row++)
	  Mj0.set_entry(row, col, M_SQRT1_2 * *it);
      }

    cout << "Mj0=" << endl << Mj0 << endl;
  }
  
  template <int d, int dT>
  void
  DSBasis<d,dT>::setup_Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp) {
    // IGPMlib reference: I_Basis_Bspline_s::Mj0ts()

    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);
    Mj0Tp.resize(njp, nj);
    
    for (unsigned int i = 0; i < MLTp.row_dimension(); i++)
      for (unsigned int k = 0; k < MLTp.column_dimension(); k++)
 	Mj0Tp.set_entry(i, k, MLTp.get_entry(i, k));
    
    for (unsigned int i = 0; i < MRTp.row_dimension(); i++)
      for (unsigned int k = 0; k < MRTp.column_dimension(); k++)
 	Mj0Tp.set_entry(njp-i-1, nj-k-1, MRTp.get_entry(i, k));
    
    int startrow = dT+ellT_l()+ell1T<d,dT>()-2*sT0;
    for (int col = dT-sT0; col < nj-(dT-sT1); col++, startrow+=2)
      {
	int row = startrow;
	for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf.aT().begin());
	     it != cdf.aT().end(); ++it, row++)
	  Mj0Tp.set_entry(row, col, M_SQRT1_2 * *it);
      }

    cout << "Mj0Tp=" << endl << Mj0Tp << endl;
  }

}

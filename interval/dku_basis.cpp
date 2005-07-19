// implementation for dku_basis.h

#include <cmath>
#include <iostream>

#include <utils/tiny_tools.h>
#include <algebra/matrix_norms.h>
#include <algebra/triangular_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/matrix_decomp.h>

#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/dm_mask.h>
#include <Rd/refinable.h>
#include <Rd/r_index.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <int d, int dT>
  typename DKUBasis<d, dT>::Index
  DKUBasis<d, dT>::firstGenerator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaLmin(), this);
  }

  template <int d, int dT>
  typename DKUBasis<d, dT>::Index
  DKUBasis<d, dT>::lastGenerator(const int j) const
  {
    assert(j >= j0());
    return Index(j, 0, DeltaRmax(j0()), this);
  }

  template <int d, int dT>
  typename DKUBasis<d, dT>::Index
  DKUBasis<d, dT>::firstWavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamin(), this);
  }

  template <int d, int dT>
  typename DKUBasis<d, dT>::Index
  DKUBasis<d, dT>::lastWavelet(const int j) const
  {
    assert(j >= j0());
    return Index(j, 1, Nablamax(j), this);
  }

  template <int d, int dT>
  DKUBasis<d, dT>::DKUBasis(DKUBiorthogonalizationMethod bio)
    : ell1_(ell1<d>()),
      ell2_(ell2<d>()),
      ell1T_(ell1T<d,dT>()),
      ell2T_(ell2T<d,dT>()),
      bio_(bio)
  {
    ellT_ = ell2T<d,dT>(); // (3.2.10)
    ell_  = ellT_-(dT-d);  // (3.2.16)

    setup_Alpha();
    setup_AlphaT();
    setup_BetaL();
    setup_BetaLT();
    setup_BetaR();
    setup_BetaRT();
    setup_GammaL();
    setup_CX_CXT();
    setup_CXA_CXAT();

    // IGPMlib reference: I_Basis_Bspline_s::Setup()

    setup_Cj();
    Matrix<double> ml = ML(); // (3.5.2)
    Matrix<double> mr = MR(); // (3.5.2)
    SparseMatrix<double> mj0;   Mj0  (ml,   mr,   mj0);   // (3.5.1)

    Matrix<double> mltp = MLTp(); // (3.5.6)
    Matrix<double> mrtp = MRTp(); // (3.5.6)
    SparseMatrix<double> mj0tp; Mj0Tp(mltp, mrtp, mj0tp); // (3.5.5)

    // construction of the wavelet basis: initial stable completion, [DKU section 4.1]
    SparseMatrix<double> FF; F(FF);         // (4.1.14)
    SparseMatrix<double> PP; P(ml, mr, PP); // (4.1.22)
    SparseMatrix<double> A, H, Hinv; GSetup(A, H, Hinv); // (4.1.1), (4.1.13)
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_Alpha() {
    // setup Alpha = (Alpha(m, r))

    const int mLow = 1-ell2T_;         // start index in (3.2.26)
    const int mUp  = 2*ellT_+ell1T_-1; // end index in (3.2.41)
    const int rUp  = dT-1;

    Alpha_.resize(mUp-mLow+1, rUp+1);

    for (int m = mLow; m <= mUp; m++)
      Alpha_(-mLow+m, 0) = 1.0; // (5.1.1)

    for (int r = 1; r <= rUp; r++) {
      double dummy = 0;
      for (int k(ell1_); k <= ell2_; k++) {
	double dummy1 = 0;
	for (int s = 0; s <= r-1; s++)
	  dummy1 += binomial(r, s) * intpower(k, r-s) * Alpha_(-mLow, s);
	dummy += cdf_.a().get_coefficient(MultiIndex<int, 1>(k)) * dummy1; // (5.1.3)
      }
      Alpha_(-mLow, r) = dummy / (ldexp(1.0, r+1) - 2.0);
    }

    for (int r(1); r <= rUp; r++) {
      for (int m = mLow; m < 0; m++) {
	double dummy = 0;
	for (int i = 0; i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * Alpha_(-mLow, r-i); // (5.1.2)
	Alpha_(-mLow+m, r) = dummy;
      }

      for (int m = 1; m <= mUp; m++) {
	double dummy = 0;
	for (int i = 0; i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * Alpha_(-mLow, r-i); // (5.1.2)
	Alpha_(-mLow+m, r) = dummy;
      }
    }
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_AlphaT() {
    // setup AlphaT

    const int mLow = 1-ell2_;        // start index in (3.2.25)
    const int mUp  = 2*ell_+ell1_-1; // end index in (3.2.40)
    const int rUp  = d-1;

    AlphaT_.resize(mUp-mLow+1, rUp+1);

    for (int m = mLow; m <= mUp; m++)
      AlphaT_(-mLow+m, 0) = 1.0; // (5.1.1)

    for (int r = 1; r <= rUp; r++) {
      double dummy = 0;
      for (int k(ell1T_); k <= ell2T_; k++) {
	double dummy1 = 0;
	for (int s = 0; s <= r-1; s++)
	  dummy1 += binomial(r, s) * intpower(k, r-s) * AlphaT_(-mLow, s); // (5.1.3)
	dummy += cdf_.aT().get_coefficient(MultiIndex<int, 1>(k)) * dummy1;
      }
      AlphaT_(-mLow, r) = dummy / (ldexp(1.0, r+1) - 2.0);
    }

    for (int r = 1; r <= rUp; r++) {
      for (int m = mLow; m < 0; m++) {
	double dummy = 0;
	for (int i = 0; i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * AlphaT_(-mLow, r-i); // (5.1.2)
	AlphaT_(-mLow+m, r) = dummy;
      }

      for (int m = 1; m <= mUp; m++) {
	double dummy = 0;
	for (int i = 0; i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * AlphaT_(-mLow, r-i); // (5.1.2)
	AlphaT_(-mLow+m, r) = dummy;
      }
    } 
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_BetaL() {
    // setup BetaL

    const int mLow = 2*ellT_+ell1T_;   // start index in (3.2.41)
    const int mUp  = 2*ellT_+ell2T_-2; // end index in (3.2.41)
    const int rUp  = dT-1;

    const int AlphamLow = 1-ell2T_;

    BetaL_.resize(mUp-mLow+1, rUp+1);

    for (int r = 0; r <= rUp; r++)
      for (int m = mLow; m <= mUp; m++) {
	double help = 0;
	for (int q = (int)ceil((m-ell2T_)/2.0); q < ellT_; q++)
	  help += Alpha_(-AlphamLow+q, r) * cdf_.aT().get_coefficient(MultiIndex<int, 1>(m-2*q)); // (3.2.31)

	BetaL_(-mLow+m, r) = help * M_SQRT1_2;
      }
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_BetaLT() {
    // setup BetaLT

    const int mLow = 2*ell_+ell1_;   // start index in (3.2.40)
    const int mUp  = 2*ell_+ell2_-2; // end index in (3.2.40)
    const int rUp  = d-1;

    const int AlphaTmLow = 1-ell2_;

    BetaLT_.resize(mUp-mLow+1, rUp+1);

    for (int r = 0; r <= rUp; r++)
      for (int m = mLow; m <= mUp; m++) {
	double help = 0;
	for (int q = (int)ceil((m-ell2_)/2.0); q < ell_; q++)
	  help += AlphaT_(-AlphaTmLow+q, r) * cdf_.a().get_coefficient(MultiIndex<int, 1>(m-2*q)); // (3.2.31)
	
	BetaLT_(-mLow+m, r) = help * M_SQRT1_2;
      }
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_BetaR() {
    // setup BetaR

    const int mLow = 2*ellT_+ell1T_;  // start index in (3.2.41)
    const int mUp  = 2*ellT_+ell2T_-2; // end index in (3.2.41)
    const int rUp  = dT-1;

    const int AlphamLow = 1-ell2T_;

    BetaR_.resize(mUp-mLow+1, rUp+1);

    for (int r = 0; r <= rUp; r++)
      for (int m = mLow; m <= mUp; m++) {
	double help = 0;
	for (int q = (int)ceil((m-ell2T_)/2.0); q < ellT_; q++)
	  help += Alpha_(-AlphamLow+q, r) * cdf_.aT().get_coefficient(MultiIndex<int, 1>(m-2*q)); // (3.2.31)

	BetaR_(-mLow+m, r) = help * M_SQRT1_2;
      }
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_BetaRT() {
    // setup BetaRT

    const int mLow = 2*ell_+ell1_;   // start index in (3.2.40)
    const int mUp  = 2*ell_+ell2_-2; // end index in (3.2.40)
    const int rUp  = d-1;

    const int AlphaTmLow = 1-ell2_;

    BetaRT_.resize(mUp-mLow+1, rUp+1);

    for (int r = 0; r <= rUp; r++)
      for (int m = mLow; m <= mUp; m++) {
	double help = 0;
	for (int q = (int)ceil((m-ell2_)/2.0); q < ell_; q++)
	  help += AlphaT_(-AlphaTmLow+q, r) * cdf_.a().get_coefficient(MultiIndex<int, 1>(m-2*q)); // (3.2.31)
	
	BetaRT_(-mLow+m, r) = help * M_SQRT1_2;
      }
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_GammaL() {
    // IGPMlib reference: imask_bspline.cpp, experimental version

//     const int llowT = ellT-dT;
//     const int lupT  = ellT-1;

    GammaL_.resize(dT, dT); // lupT-llowT+1 each
	
    // 1. compute the integrals
    //      z(s,t) = \int_0^1\phi(x-s)\tilde\phi(x-t)\,dx
    //    exactly with the [DM] trick
    MultivariateRefinableFunction<DMMask2<HaarMask, CDFMask_primal<d>, CDFMask_dual<d, dT> >, 2> zmask;
    InfiniteVector<double, MultiIndex<int, 2> > zvalues(zmask.evaluate());

    // 2. compute the integrals
    //      I(nu,mu) = \int_0^\infty\phi(x-\nu)\tilde\phi(x-\mu)\,dx
    //    exactly using the z(s,t) values

    const int I1Low = -ell2_+1;
    const int I1Up  = ell_-d+dT-1;
    const int I2Low = -ell2T_+1;
    const int I2Up  = ellT_-1;

    Matrix<double> I(I1Up-I1Low+1, I2Up-I2Low+1);
    for (int nu = I1Low; nu < -ell1_; nu++)
      for (int mu = I2Low; mu < -ell1T_; mu++) {
	double help(0);
	int diff(mu - nu);
	for (int s(-ell2_+1); s <= nu; s++)
	  help += zvalues.get_coefficient(MultiIndex<int, 2>(s, s+diff));
	I(-I1Low+nu, -I2Low+mu) = help; // (5.1.7)
      }
    for (int nu = -I1Low; nu <= I1Up; nu++)
      for (int mu = I2Low; mu <= I2Up; mu++) {
	if ((nu >= -ell1_) || ((nu <= ell_-1) && (mu >= -ell1T_)))
	  I(-I1Low+nu, -I2Low+mu) = (nu == mu ? 1 : 0); // (5.1.6)
      }

    // 3. finally, compute the Gramian GammaL
    //    (offsets: entry (r,k) <-> index (r+ellT()-dT,k+ellT()-dT) )

    const int AlphamLow = 1-ell2T_;
    const int AlphaTmLow = 1-ell2_;

    for (int r(0); r < d; r++)
      for (int k(0); k < dT; k++) {
	double help(0);
	for (int nu = I1Low; nu < ell_; nu++)
	  for (int mu = I2Low; mu <= I2Up; mu++)
	    help += AlphaT_(-AlphaTmLow+nu, r) * Alpha_(-AlphamLow+mu, k) * I(-I1Low+nu, -I2Low+mu);
	GammaL_(r, k) = help; // (5.1.4)
      }
    for (int r(d); r <= dT-1; r++)
      for (int k(0); k <= dT-1; k++) {
	double help(0);
	for (int mu = I2Low; mu <= I2Up; mu++)
	  help += Alpha_(-AlphamLow+mu, k) * I(-I1Low+ell_-d+r, -I2Low+mu);
	GammaL_(r, k) = help; // (5.1.5)
      }
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_CX_CXT() {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

//     const int llowT = ellT_-dT;
//     const int lupT  = ellT_-1;

    CL_.resize(dT, dT);  // lupT-llowT+1 each, same offsets as GammaL
    CLT_.resize(dT, dT); // "
    CR_.resize(dT, dT);  // "
    CRT_.resize(dT, dT); // "

    // setup CL and CLT
    // (offsets: entry (i,j) <-> index (i+ellT()-dT,j+ellT()-dT) )
    //
    if (bio_ == none) {
      for (unsigned int i(0); i < dT; i++)
	CL_(i, i) = 1.0; // identity matrix
      
      Matrix<double> CLGammaLInv;
      QRDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
      CLT_ = transpose(CLGammaLInv);
    }
    
    if (bio_ == SVD) {
      MathTL::SVD<double> svd(GammaL_);
      Matrix<double> U, V;
      Vector<double> S;
      svd.getU(U);
      svd.getV(V);
      svd.getS(S);
      
      for (unsigned int i(0); i < dT; i++) {
	S[i] = 1.0 / sqrt(S[i]);
	for (unsigned int j(0); j < dT; j++) {
	  CL_(i, j) = S[i] * U(j, i);
	  CLT_(i, j) = S[i] * V(i, j);
	}
      }

      CL_.compress(1e-12);
      CLT_.compress(1e-12);
    }
    
    if (bio_ == Bernstein) {
      double b(0);
      if (d == 2)
	b = 0.7; // cf. [DKU]
      else {
	b = 3.6;
      }
      
      // CL_ : transformation to Bernstein polynomials (triangular)
      //   (Z_b^m)_{r,l} = (-1)^{r-1}\binom{m-1}{r}\binom{r}{l}b^{-r}, r>=l, 0 otherwise
      for (unsigned int j(0); j < d; j++)
	for (unsigned int k(0); k <= j; k++)
	  CL_(k, j) = minus1power(j-k) * std::pow(b, -(int)j) * binomial(d-1, j) * binomial(j, k);
      for (unsigned int j(d); j < dT; j++)
	CL_(j, j) = 1.0;
      
      Matrix<double> CLGammaLInv;
      QRDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
      CLT_ = transpose(CLGammaLInv);
      
      CLT_.compress(1e-12);
    }
    
    if (bio_ == partialSVD) {
      MathTL::SVD<double> svd(GammaL_);
      Matrix<double> U, V;
      Vector<double> S;
      svd.getU(U);
      svd.getV(V);
      svd.getS(S);
      
      Matrix<double> GammaLInv;
      QRDecomposition<double>(GammaL_).inverse(GammaLInv);
      const double a = 1.0 / GammaLInv(0, 0);
      
      Matrix<double> R(dT, dT);
      for (unsigned int i(0); i < dT; i++)
	S[i] = sqrt(S[i]);
      for (unsigned int j(0); j < dT; j++)
	R(0, j) = a * V(0, j) / S[j];
      for (unsigned int i(1); i < dT; i++)
	for (unsigned int j(0); j < dT; j++)
	  R(i, j) = U(j, i) * S[j];
      
      for (unsigned int i(0); i < dT; i++)
	for (unsigned int j(0); j < dT; j++)
	  U(i, j) /= S[i];
      CL_ = R*U;

      CL_.compress(1e-12);

      Matrix<double> CLGammaLInv;
      QRDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
      CLT_ = transpose(CLGammaLInv);

      CLT_.compress(1e-12);
    }

    if (bio_ == BernsteinSVD) {
      double b(0);
      if (d == 2)
	b = 0.7; // cf. [DKU]
      else {
	b = 3.6;
      }
      
      // CL_ : transformation to Bernstein polynomials (triangular)
      //   (Z_b^m)_{r,l} = (-1)^{r-1}\binom{m-1}{r}\binom{r}{l}b^{-r}, r>=l, 0 otherwise
      for (unsigned int j(0); j < d; j++)
	for (unsigned int k(0); k <= j; k++)
	  CL_(k, j) = minus1power(j-k) * std::pow(b, -(int)j) * binomial(d-1, j) * binomial(j, k);
      for (unsigned int j(d); j < dT; j++)
	CL_(j, j) = 1.0;
      
      Matrix<double> GammaLNew(CL_ * GammaL_);
      
      MathTL::SVD<double> svd(GammaLNew);
      Matrix<double> U, V;
      Vector<double> S;
      svd.getU(U);
      svd.getV(V);
      svd.getS(S);
      
      Matrix<double> GammaLNewInv;
      QRDecomposition<double>(GammaLNew).inverse(GammaLNewInv);
      const double a = 1.0 / GammaLNewInv(0, 0);
      
      Matrix<double> R(dT, dT);
      for (unsigned int i(0); i < dT; i++)
	S[i] = sqrt(S[i]);
      for (unsigned int j(0); j < dT; j++)
	R(0, j) = a * V(0, j) / S[j];
      for (unsigned int i(1); i < dT; i++)
	for (unsigned int j(0); j < dT; j++)
	  R(i, j) = U(j, i) * S[j];
      
      for (unsigned int i(0); i < dT; i++)
	for (unsigned int j(0); j < dT; j++)
	  U(i, j) /= S[i];
      CL_ = R*U*CL_;
      
      CL_.compress(1e-12);

      Matrix<double> CLGammaLInv;
      QRDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
      CLT_ = transpose(CLGammaLInv);

      CLT_.compress(1e-12);
    }
    
#if 0
    // check biorthogonality of the matrix product CL * GammaL * (CLT)^T
    Matrix<double> check(CL_ * GammaL_ * transpose(CLT_));
    for (unsigned int i(0); i < check.row_dimension(); i++)
      check(i, i) -= 1;
    cout << "error for CLT: " << row_sum_norm(check) << endl;
#endif

    QRDecomposition<double>(CL_).inverse(inv_CL_);
    QRDecomposition<double>(CLT_).inverse(inv_CLT_);

    CL_.mirror(CR_);
    CLT_.mirror(CRT_);

    QRDecomposition<double>(CR_).inverse(inv_CR_);
    QRDecomposition<double>(CRT_).inverse(inv_CRT_);
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_CXA_CXAT() {
    // IGPMlib reference: I_Mask_Bspline::EvalCL(), ::EvalCR()

    const int llklow  = 1-ell2_;   // offset 1 in CLA (and AlphaT), see (3.2.25)
    const int llkup   = ell_-1;
    const int llklowT = 1-ell2T_;  // offset 1 in CLAT (and Alpha), see (3.2.26)
    const int llkupT  = ellT_-1;

    const int llow    = ell_-d;    // offset 1 in CL (and CLT), offset 2 in AlphaT
    const int lup     = ell_-1;

    const int llowT   = ellT_-dT;  // == llow, offset 2 in CLA and CLAT
    const int lupT    = ellT_-1;

    // setup CLA <-> AlphaT * (CL)^T
    CLA_.resize(llkup-llklow+1, lupT-llowT+1);

    for (int i = llklow; i <= llkup; i++) // the (3.2.25) bounds
      for (int r = llowT; r <= lupT; r++) {
	double help = 0;
	for (int m = llow; m <= lup; m++)
	  help += CL_(-llowT+r, -llowT+m) * AlphaT_(-llklow+i, -llow+m);
	CLA_(-llklow+i, -llowT+r) = help;
      }
    for (int i = lup+1; i <= llkup; i++)
      for (int r = llowT; r <= lupT; r++)
	CLA_(-llklow+i, -llowT+r) += CL_(-llowT+r, -llowT+i);

    CLA_.compress(1e-12);

    CLA_.mirror(CRA_);

    // setup CLAT <-> Alpha * (CLT)^T
    CLAT_.resize(llkupT-llklowT+1, lupT-llowT+1);
    for (int i = llklowT; i <= llkupT; i++) // the (3.2.26) bounds
      for (int r = llowT; r <= lupT; r++) {
	double help = 0;
	for (int m = llowT; m <= lupT; m++)
	  help += CLT_(-llowT+r, -llowT+m) * Alpha_(-llklowT+i, -llowT+m);
  	CLAT_(-llklowT+i, -llowT+r) = help;
      }

    CLAT_.compress(1e-12);

    CLAT_.mirror(CRAT_);
  }

  template <int d, int dT>
  void DKUBasis<d, dT>::setup_Cj() {
    // IGPMlib reference: I_Basis_Bspline_s::setup_Cj(), ::put_Mat()

    // (5.2.5)
    Cj_.diagonal(Deltasize(j0()), 1.0);
    Cj_.set_block(0, 0, CL_);
    Cj_.set_block(Deltasize(j0())-dT, Deltasize(j0())-dT, CR_);

    inv_Cj_.diagonal(Deltasize(j0()), 1.0);
    inv_Cj_.set_block(0, 0, inv_CL_);
    inv_Cj_.set_block(Deltasize(j0())-dT, Deltasize(j0())-dT, inv_CR_);

    CjT_.diagonal(Deltasize(j0()), 1.0);
    CjT_.set_block(0, 0, CLT_);
    CjT_.set_block(Deltasize(j0())-dT, Deltasize(j0())-dT, CRT_);

    inv_CjT_.diagonal(Deltasize(j0()), 1.0);
    inv_CjT_.set_block(0, 0, inv_CLT_);
    inv_CjT_.set_block(Deltasize(j0())-dT, Deltasize(j0())-dT, inv_CRT_);

    Cjp_.diagonal(Deltasize(j0()+1), 1.0);
    Cjp_.set_block(0, 0, CL_);
    Cjp_.set_block(Deltasize(j0()+1)-dT, Deltasize(j0()+1)-dT, CR_);

    inv_Cjp_.diagonal(Deltasize(j0()+1), 1.0);
    inv_Cjp_.set_block(0, 0, inv_CL_);
    inv_Cjp_.set_block(Deltasize(j0()+1)-dT, Deltasize(j0()+1)-dT, inv_CR_);

    CjpT_.diagonal(Deltasize(j0()+1), 1.0);
    CjpT_.set_block(0, 0, CLT_);
    CjpT_.set_block(Deltasize(j0()+1)-dT, Deltasize(j0()+1)-dT, CRT_);

    inv_CjpT_.diagonal(Deltasize(j0()+1), 1.0);
    inv_CjpT_.set_block(0, 0, inv_CLT_);
    inv_CjpT_.set_block(Deltasize(j0()+1)-dT, Deltasize(j0()+1)-dT, inv_CRT_);  
  }
  
  template <int d, int dT>
  Matrix<double>
  DKUBasis<d, dT>::ML() const
  {
    // IGPMlib reference: I_Basis_Bspline_s::ML()

    const int llow = ell_-d;           // the (3.5.2) bounds
    const int lup  = ell_-1;
    const int MLup = ell2_+2*ell_-2;
    const int AlphaToffset = 1-ell2_;
    const int BetaLToffset = 2*ell_+ell1_;

    Matrix<double> ML(MLup-llow+1, lup-llow+1);

    for (int row = 0; row < d; row++)
      ML(row, row) = 1.0 / sqrt(ldexp(1.0, 2*row+1));

    for (int m = ell_; m <= 2*ell_+ell1_-1; m++)
      for (int k = 0; k < d; k++)
    	ML(-llow+m, k) = AlphaT_(-AlphaToffset+m, k) / sqrt(ldexp(1.0, 2*k+1));
    
    for (int m = 2*ell_+ell1_; m <= MLup; m++)
      for (int k = 0; k < d; k++)
  	ML(-llow+m, k) = BetaLT_(-BetaLToffset+m, k);

    return ML;
  }

  template <int d, int dT>
  Matrix<double>
  DKUBasis<d, dT>::MR() const
  {
    // IGPMlib reference: I_Basis_Bspline_s::MR()
    // (remark: in this form identical to ML(), in I_Basis_Bspline_s they
    // use elll and ellr, which is identical here)

    const int llow = ell_-d;           // the (3.5.2) bounds
    const int lup  = ell_-1;
    const int MRup  = ell2_+2*ell_-2;
    const int AlphaToffset = 1-ell2_;
    const int BetaRToffset = 2*ell_+ell1_;

    Matrix<double> MR(MRup-llow+1, lup-llow+1);

    for (int row = 0; row < d; row++)
      MR(row, row) = 1.0 / sqrt(ldexp(1.0, 2*row+1));

    for (int m = ell_; m <= 2*ell_+ell1_-1; m++)
      for (int k = 0; k < d; k++)
	MR(-llow+m, k) = AlphaT_(-AlphaToffset+m, k) / sqrt(ldexp(1.0, 2*k+1));
    
    for (int m = 2*ell_+ell1_; m <= MRup; m++)
      for (int k = 0; k < d; k++)
  	MR(-llow+m, k) = BetaRT_(-BetaRToffset+m, k);

    return MR;
  }

  template <int d, int dT>
  Matrix<double>
  DKUBasis<d, dT>::MLTp() const
  {
    // IGPMlib reference: I_Basis_Bspline_s::MLts()

    const int llowT = ellT_-dT;           // the (3.5.6) bounds
    const int lupT  = ellT_-1;
    const int MLTup = ell2T_+2*ellT_-2;
    const int Alphaoffset = 1-ell2T_;
    const int BetaLoffset = 2*ellT_+ell1T_;

    Matrix<double> MLTp(MLTup-llowT+1, lupT-llowT+1);

    for (int row = 0; row < dT; row++)
      MLTp(row, row) = 1.0 / sqrt(ldexp(1.0, 2*row+1));

    for (int m = ellT_; m <= 2*ellT_+ell1T_-1; m++)
      for (int k = 0; k < d; k++)
    	MLTp(-llowT+m, k) = Alpha_(-Alphaoffset+m, k) / sqrt(ldexp(1.0, 2*k+1));
    
    for (int m = 2*ellT_+ell1T_; m <= MLTup; m++)
      for (int k = 0; k < d; k++)
  	MLTp(-llowT+m, k) = BetaL_(-BetaLoffset+m, k);

    return MLTp;
  }

  template <int d, int dT>
  Matrix<double>
  DKUBasis<d, dT>::MRTp() const
  {
    // IGPMlib reference: I_Basis_Bspline_s::MRts()

    const int llowT = ellT_-dT;           // the (3.5.6) bounds
    const int lupT  = ellT_-1;
    const int MRTup = ell2T_+2*ellT_-2;
    const int Alphaoffset = 1-ell2T_;
    const int BetaLoffset = 2*ellT_+ell1T_;

    Matrix<double> MRTp(MRTup-llowT+1, lupT-llowT+1);

    for (int row = 0; row < dT; row++)
      MRTp(row, row) = 1.0 / sqrt(ldexp(1.0, 2*row+1));

    for (int m = ellT_; m <= 2*ellT_+ell1T_-1; m++)
      for (int k = 0; k < d; k++)
    	MRTp(-llowT+m, k) = Alpha_(-Alphaoffset+m, k) / sqrt(ldexp(1.0, 2*k+1));
    
    for (int m = 2*ellT_+ell1T_; m <= MRTup; m++)
      for (int k = 0; k < d; k++)
  	MRTp(-llowT+m, k) = BetaL_(-BetaLoffset+m, k);

    return MRTp;
  }

  template <int d, int dT>
  void
  DKUBasis<d, dT>::Mj0(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& Mj0)
  {
    // IGPMlib reference: I_Basis_Bspline_s::Mj0()
    
    // TODO: enhance readability! (<-> [DKU section 3.5])

    int p = (1 << j0()) - 2*ell_ - (d%2) + 1;
    int q = (1 << j0()) - 4*ell_ - 2*(d%2) + d + 1; // this value is overwritten below

    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);

    Mj0.resize(njp, nj);

    const int alowc = d+1;
    const int aupc  = d+p;
    const int alowr = d+ell_+ell1_;
    
    for (int i = 0; i < (int)ML.row_dimension(); i++)
      for (int k = 0; k < (int)ML.column_dimension(); k++)
	Mj0.set_entry(i, k, ML.get_entry(i, k));

    for (int i = 0; i < (int)MR.row_dimension(); i++)
      for (int k = 0; k < (int)MR.column_dimension(); k++)
	Mj0.set_entry(njp-i-1, nj-k-1, MR.get_entry(i, k));

    q = alowr;
    for (int r = alowc; r <= aupc; r++)
      {
	p = q;
	for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf_.a().begin());
	     it != cdf_.a().end(); ++it)
	  {
	    p++;
	    Mj0.set_entry(p-1, r-1, M_SQRT1_2 * *it); // quick hack, TODO: eliminate the "-1"
	  }
	q += 2;
      }
    
    // after this routine: offsets in Mj0 are both llow == ell-d
  }

  template <int d, int dT>
  void
  DKUBasis<d, dT>::Mj0Tp(const Matrix<double>& MLTp, const Matrix<double>& MRTp, SparseMatrix<double>& Mj0Tp)
  {
    // IGPMlib reference: I_Basis_Bspline_s::Mj0ts()

    // TODO: enhance readability! (<-> [DKU section 3.5])
    
    int p = (1 << j0()) - 2*ellT_ - (dT%2) + 1;
    int q = (1 << j0()) - 4*ellT_ - 2*(dT%2) + dT + 1; // this value is overwritten below

    const int nj  = Deltasize(j0());
    const int njp = Deltasize(j0()+1);

    Mj0Tp.resize(njp, nj);

    const int atlowc = dT+1;
    const int atupc  = dT+p;
    const int atlowr = dT+ellT_+ell1T_;
    
    for (int i = 0; i < (int)MLTp.row_dimension(); i++)
      for (int k = 0; k < (int)MLTp.column_dimension(); k++)
	Mj0Tp.set_entry(i, k, MLTp.get_entry(i, k));

    for (int i = 0; i < (int)MRTp.row_dimension(); i++)
      for (int k = 0; k < (int)MRTp.column_dimension(); k++)
	Mj0Tp.set_entry(njp-i-1, nj-k-1, MRTp.get_entry(i, k));

    q = atlowr;
    for (int r = atlowc; r <= atupc; r++)
      {
	p = q+1;
	for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf_.aT().begin());
	     it != cdf_.aT().end(); ++it, p++)
	  {
	    Mj0Tp.set_entry(p-1, r-1, M_SQRT1_2 * *it); // quick hack, TODO: eliminate the "-1"
	  }
	q += 2;
      }
    
    // after this routine: offsets in Mj0Tp are both llowT == ellT-dT
  }

  template <int d, int dT>
  void
  DKUBasis<d, dT>::F(SparseMatrix<double>& FF)
  {
    // IGPMlib reference: I_Basis_Bspline_s::F()

    const int FLow = ell_+(d%2);       // start column index for F_j in (4.1.14)
    const int FUp  = (1<<j0()) - ell_; // end column index for F_j in (4.1.14)

    // (4.1.14):
    FF.resize(Deltasize(j0()+1), 1<<j0());
    for (int r = 1; r <= FLow-1; r++)
      FF.set_entry(r+d-1, r-1, 1.0);
    
    int i = d+ell_+(d%2)-1;
    for (int r = FLow; r <= FUp; r++)
      {
	FF.set_entry(i, r-1, 1.0);
	i += 2;
      } 

    i = Deltasize(j0()+1)-d-1;
    for (int r = 1<<j0(); r >= FUp+1; r--)
      {
	FF.set_entry(i, r-1, 1.0);
	i--;
      }
  }

  template <int d, int dT>
  void
  DKUBasis<d, dT>::P(const Matrix<double>& ML, const Matrix<double>& MR, SparseMatrix<double>& PP)
  {
    // IGPMlib reference: I_Basis_Bspline_s::P()
    
    // (4.1.22):
    PP.diagonal(Deltasize(j0()+1), 1.0);
    
    for (int i = 0; i < (int)ML.row_dimension(); i++)
      for (int k = 0; k < (int)ML.column_dimension(); k++)
	PP.set_entry(i, k, ML.get_entry(i, k));

    for (int i = 0; i < (int)MR.row_dimension(); i++)
      for (int k = 0; k < (int)MR.column_dimension(); k++)
	PP.set_entry(Deltasize(j0()+1)-i-1, Deltasize(j0()+1)-k-1, MR.get_entry(i, k));
  }

  template <int d, int dT>
  void
  DKUBasis<d, dT>::GSetup(SparseMatrix<double>& A, SparseMatrix<double>& H, SparseMatrix<double>& Hinv)
  {
    // IGPMlib reference: I_Basis_Bspline_s::GSetup()

    // (4.1.13):
    A.resize(Deltasize(j0()+1), Deltasize(j0()));

    for (int r = 0; r < d; r++)
      A.set_entry(r, r, 1.0);

    for (int r = d, q = d+ell_+ell1_; r < Deltasize(j0())-d; r++)
      {
 	int p = q;

 	for (MultivariateLaurentPolynomial<double, 1>::const_iterator it(cdf_.a().begin());
 	     it != cdf_.a().end(); ++it)
 	  {
 	    A.set_entry(p, r, M_SQRT1_2 * *it);
	    p++;
 	  }

	q += 2;
      }

    for (int r = Deltasize(j0()+1)-d, q = Deltasize(j0())-d; r < Deltasize(j0()+1); r++, q++)
      A.set_entry(r, q, 1.0);

//     H.Redimension(Deltajp1, Deltajp1);
//     H.Identity(1.0);

//     inverseH.Redimension(Deltajp1, Deltajp1);
//     inverseH.Identity(1.0);

//     H.setindexr(_ll-_d);
//     H.setindexc(_ll-_d);

//     inverseH.setindexr(_ll-_d);
//     inverseH.setindexc(_ll-_d);

  }

  template <int d, int dT>
  SampledMapping<1>
  DKUBasis<d, dT>::evaluate(const Index& lambda,
			    const bool primal,
			    const int resolution) const
  {
    if (lambda.e() == 0) { // generator
      if (primal) {
	if (lambda.k() <= DeltaLTmax()) {
	    // left boundary generator
	    InfiniteVector<double, RIndex> coeffs;
	    for(int i(0); i < ellT_+ell2_-1; i++) {
	      double v(CLA_(i, lambda.k()-ellT_+dT));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, i+1-ell2_), v);
	    }
	    return cdf_.evaluate(0, coeffs, primal, 0, 1, resolution);
	} else {
	  if (lambda.k() >= DeltaRTmin(lambda.j())) {
	    // right boundary generator
	    InfiniteVector<double, RIndex> coeffs;
	    for (int i(0); i < ellT_+ell2_-1; i++) {
	      double v(CRA_(i, lambda.k()+dT-1-DeltaRmax(lambda.j())));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, DeltaRmax(lambda.j())-ellT_-ell2_+2+i), v);
	    }
	    return cdf_.evaluate(0, coeffs, primal, 0, 1, resolution);
	  } else {
	    // inner generator
	    return cdf_.evaluate(0, RIndex(lambda.j(), 0, lambda.k()),
				 primal, 0, 1, resolution);
	  }
	}
      } else {
	// dual
	if (lambda.k() <= DeltaLTmax()) {
	  // left boundary generator
	  InfiniteVector<double, RIndex> coeffs;
	  for (int i(0); i < ellT_+ell2T_-1; i++) {
	    double v(CLAT_(i, lambda.k()-ell_+d));
	    if (v != 0)
	      coeffs.set_coefficient(RIndex(lambda.j(), 0, i+1-ell2T_), v);
	  }
	  return cdf_.evaluate(0, coeffs, primal, 0, 1, resolution);
	} else {
	  if (lambda.k() >= DeltaRTmin(lambda.j())) {
	    // right boundary generator
	    InfiniteVector<double, RIndex> coeffs;
	    for (int i(0); i < ellT_+ell2T_-1; i++) {
	      double v(CRAT_(i, lambda.k()+dT-1-DeltaRmax(lambda.j())));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, DeltaRmax(lambda.j())-ellT_-ell2_+2+i), v);
	    }
	    return cdf_.evaluate(0, coeffs, primal, 0, 1, resolution);
	  } else {
	    // inner generator
	    return cdf_.evaluate(0, RIndex(lambda.j(), 0, lambda.k()),
				 primal, 0, 1, resolution);
	  }
	}
      }
    } else {
      // expand wavelet as generators of a higher scale
      // TODO!!!
    }
    
    return SampledMapping<1>(); // dummy return for the compiler
  }

//   template <int d, int dT>
//   SampledMapping<1> evaluate(const InfiniteVector<double, Index>& coeffs,
// 			     const bool primal,
// 			     const int resolution) const;


}

// implementation for dku.h

#include <cmath>
#include <iostream>
#include <utils/tiny_tools.h>
#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/dm_mask.h>
#include <Rd/multi_refinable.h>
#include <algebra/matrix_norms.h>
#include <algebra/triangular_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/matrix_decomp.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  // helper routine to check Alpha
  template <int d, int dT>
  double DKUBasis<d, dT>::alpha_(const int m, const int r) const
  {
    double result(0);

    if (r == 0)
      result = 1.0; // (5.1.1)
    else
      {
	if (m == 0)
	  {
	    double dummy(0);
	    for (int k(ell1_); k <= ell2_; k++)
	      {
		double dummy1(0);
		for (int s(0); s <= r-1; s++)
		  dummy1 += binomial(r, s) * intpower(k, r-s) * alpha_(0, s);
		dummy += cdf_.a().get_coefficient(k) * dummy1; // (5.1.3)
	      }
	    result = dummy / (ldexp(1.0, r+1) - 2.0);
	  }
	else
	  {
  	    for (int i(0); i <= r; i++)
 	      result += binomial(r, i) * intpower(m, i) * alpha_(0, r-i); // (5.1.2)
	  }
      }

    return result;
  }

  // helper routine to check AlphaT
  template <int d, int dT>
  double DKUBasis<d, dT>::alphaT_(const int m, const int r) const
  {
    double result(0);

    if (r == 0)
      result = 1.0; // (5.1.1)
    else
      {
	if (m == 0)
	  {
	    double dummy(0);
	    for (int k(ell1T_); k <= ell2T_; k++)
	      {
		double dummy1(0);
		for (int s(0); s <= r-1; s++)
		  dummy1 += binomial(r, s) * intpower(k, r-s) * alphaT_(0, s);
		dummy += cdf_.aT().get_coefficient(k) * dummy1; // (5.1.3)
	      }
	    result = dummy / (ldexp(1.0, r+1) - 2.0);
	  }
	else
	  {
 	    for (int i(0); i <= r; i++)
	      result += binomial(r, i) * intpower(m, i) * alphaT_(0, r-i); // (5.1.2)
	  }
      }

    return result;
  }

  // helper routine to check BetaL
  template <int d, int dT>
  double DKUBasis<d, dT>::betaL_(const int m, const int r) const
  {
    double result(0);

    for (int q((int)ceil((m-ell2T_)/2.0)); q <= ellT_-1; q++)
      result += Alpha_(q, r) * cdf_.aT().get_coefficient(m-2*q); // (3.2.31)

    return result * M_SQRT1_2;
  }

  // helper routine to check BetaL
  template <int d, int dT>
  double DKUBasis<d, dT>::betaLT_(const int m, const int r) const
  {
    double result(0);

    for (int q((int)ceil((m-ell2_)/2.0)); q <= ell_-1; q++)
      result += AlphaT_(q, r) * cdf_.a().get_coefficient(m-2*q); // (3.2.31)

    return result * M_SQRT1_2;
  }

  template <int d, int dT>
  DKUBasis<d, dT>::DKUBasis(DKUBiorthogonalizationMethod bio)
    : ell1_(ell1<d>()),
      ell2_(ell2<d>()),
      ell1T_(ell1T<d,dT>()),
      ell2T_(ell2T<d,dT>()),
      GammaL_(dT, dT),
      CL_(dT, dT),
      CLT_(dT, dT)
  {
    ellT_ = ell2T<d,dT>(); // (3.2.10)
    ell_ = ellT_-(dT-d); // (3.2.16)

    // setup Alpha
    Alpha_.resize(ellT_+ell2T_-1, dT);
    for (unsigned int m(0); m < Alpha_.row_dimension(); m++)
      Alpha_(m, 0) = 1.0; // (5.1.1)
    for (int r(1); r < (int)Alpha_.column_dimension(); r++) {
      double dummy(0);
      for (int k(ell1_); k <= ell2_; k++) {
	double dummy1(0);
	for (int s(0); s <= r-1; s++)
	  dummy1 += binomial(r, s) * intpower(k, r-s) * Alpha_(ell2T_-1, s);
	dummy += cdf_.a().get_coefficient(k) * dummy1; // (5.1.3)
      }
      Alpha_(ell2T_-1, r) = dummy / (ldexp(1.0, r+1) - 2.0);
    }
    for (int r(1); r < (int)Alpha_.column_dimension(); r++) {
      for (int m(-ell2T_+1); m < 0; m++) {
	double dummy(0);
	for (int i(0); i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * Alpha_(ell2T_-1, r-i); // (5.1.2)
	Alpha_(m+ell2T_-1, r) = dummy;
      }
      for (int m(1); m <= ellT_-1; m++) {
	double dummy(0);
	for (int i(0); i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * Alpha_(ell2T_-1, r-i); // (5.1.2)
	Alpha_(m+ell2T_-1, r) = dummy;
      }
    }
    
    // setup AlphaT
    AlphaT_.resize(ell_+ell2_-1, d);
    for (unsigned int m(0); m < AlphaT_.row_dimension(); m++)
      AlphaT_(m, 0) = 1.0; // (5.1.1)
    for (int r(1); r < (int)AlphaT_.column_dimension(); r++) {
      double dummy(0);
      for (int k(ell1T_); k <= ell2T_; k++) {
	double dummy1(0);
	for (int s(0); s <= r-1; s++)
	  dummy1 += binomial(r, s) * intpower(k, r-s) * AlphaT_(ell2_-1, s); // (5.1.3)
	dummy += cdf_.aT().get_coefficient(k) * dummy1;
      }
      AlphaT_(ell2_-1, r) = dummy / (ldexp(1.0, r+1) - 2.0);
    }
    for (int r(1); r < (int)AlphaT_.column_dimension(); r++) {
      for (int m(-ell2_+1); m < 0; m++) {
	double dummy(0);
	for (int i(0); i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * AlphaT_(ell2_-1, r-i); // (5.1.2)
	AlphaT_(m+ell2_-1, r) = dummy;
      }
      for (int m(1); m <= ell_-1; m++) {
	double dummy(0);
	for (int i(0); i <= r; i++)
	  dummy += binomial(r, i) * intpower(m, i) * AlphaT_(ell2_-1, r-i); // (5.1.2)
	AlphaT_(m+ell2_-1, r) = dummy;
      }
    }
    
    // setup BetaL
    BetaL_.resize(ell2T_-ell1T_-1, dT);
    for (int r(0); r < dT; r++)
      for (int m(2*ellT_+ell1T_); m <= 2*ellT_+ell2T_-2; m++)
	BetaL_(m-2*ellT_-ell1T_, r) = betaL_(m, r);

    // setup BetaLT
    BetaLT_.resize(ell2_-ell1_-1, d);
    for (int r(0); r < d; r++)
      for (int m(2*ell_+ell1_); m <= 2*ell_+ell2_-2; m++)
	BetaLT_(m-2*ell_-ell1_, r) = betaLT_(m, r);
    
    // setup GammaL:
    //
    // 1. compute the integrals
    //      z(s,t) = \int_0^1\phi(x-s)\tilde\phi(x-t)\,dx
    //    exactly with the [DM] trick
    MultivariateRefinableFunction<DMMask2<HaarMask, CDFMask_primal<d>, CDFMask_dual<d, dT> >, 2> zmask;
    InfiniteVector<double, MultiIndex<int, 2> > zvalues(zmask.evaluate());
    //
    // 2. compute the integrals
    //      I(nu,mu) = \int_0^\infty\phi(x-\nu)\tilde\phi(x-\mu)\,dx
    //    exactly using the z(s,t) values
    Matrix<double> I(ell_-d+dT-1+ell2_, ellT_-1+ell2T_);
    for (int nu(-ell2_+1); nu <= -ell1_-1; nu++)
      for (int mu(-ell2T_+1); mu <= -ell1T_-1; mu++) {
	double help(0);
	int diff(mu - nu);
	for (int s(-ell2_+1); s <= nu; s++)
	  help += zvalues.get_coefficient(MultiIndex<int, 2>(s, s+diff));
	I(nu+ell2_-1, mu+ell2T_-1) = help; // (5.1.7)
      }
    for (int nu(-ell2_+1); nu <= ell_-d+dT-1; nu++)
      for (int mu(-ell2T_+1); mu <= ellT_-1; mu++) {
	if ((nu >= -ell1_) || ((nu <= ell_-1) && (mu >= -ell1T_)))
	  I(nu+ell2_-1, mu+ell2T_-1) = (nu == mu ? 1 : 0); // (5.1.6)
      }
    //
    // 3. finally, compute the Gramian GammaL
    for (int r(0); r <= d-1; r++)
      for (int k(0); k <= dT-1; k++) {
	double help(0);
	for (int nu(-ell2_+1); nu <= ell_-1; nu++)
	  for (int mu(-ell2T_+1); mu <= ellT_-1; mu++)
	    help += AlphaT_(nu+ell2_-1, r) * Alpha_(mu+ell2T_-1, k) * I(nu+ell2_-1, mu+ell2T_-1);
	GammaL_(r, k) = help; // (5.1.4)
      }
    for (int r(d); r <= dT-1; r++)
      for (int k(0); k <= dT-1; k++) {
	double help(0);
	for (int mu(-ell2T_+1); mu <= ellT_-1; mu++)
	  help += Alpha_(mu+ell2T_-1, k) * I(ell_-d+r+ell2_-1, mu+ell2T_-1);
	GammaL_(r, k) = help; // (5.1.5)
      }
    
    // setup CL and CLT
    if (bio == none) {
      for (unsigned int i(0); i < dT; i++)
	CL_(i, i) = 1.0;

      Matrix<double> CLGammaLInv;
      QRDecomposition<double>(CL_ * GammaL_).inverse(CLGammaLInv);
      CLT_ = transpose(CLGammaLInv);
    }

    if (bio == SVD) {
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
    }

    if (bio == Bernstein) {
      const double b(1);
      const double bT(1);

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
    }

#if 0
      // setup basis transformation matrices
      //   (Z_b^m)_{r,l} = (-1)^{r-1}\binom{m-1}{r}\binom{r}{l}b^{-r}, r>=l, 0 otherwise
      LowerTriangularMatrix<double> Z(d, d), ZT(dT, dT);
      for (unsigned int j(0); j < d; j++)
	for (unsigned int k(0); k <= j; k++)
	  Z(j, k) = minus1power(j-k) * std::pow(b, -(int)j) * binomial(d-1, j) * binomial(j, k);
      for (unsigned int j(0); j < dT; j++)
	for (unsigned int k(0); k <= j; k++)
	  ZT(j, k) = minus1power(j-k) * std::pow(bT, -(int)j) * binomial(dT-1, j) * binomial(j, k);      

      // setup backtransformation matrices
      LowerTriangularMatrix<double> Zinv(d, d), ZTinv(dT, dT);
      for (unsigned int j(0); j < d; j++)
	for (unsigned int k(0); k <= j; k++)
	  Zinv(j, k) = std::pow(b, (int)k) / binomial(d-1, k) * binomial(j, k);
      for (unsigned int j(0); j < dT; j++)
	for (unsigned int k(0); k <= j; k++)
	  ZTinv(j, k) = std::pow(bT, (int)k) / binomial(dT-1, k) * binomial(j, k);
#endif



#if 1
    // check biorthogonality of matrix product CL * GammaL * (CLT)^T
    Matrix<double> check(CL_ * GammaL_ * transpose(CLT_));
    for (unsigned int i(0); i < check.row_dimension(); i++)
      check(i, i) -= 1;
    cout << "error for CLT: " << row_sum_norm(check) << endl;
#endif


  }
}

// implementation for dku.h

#include <cmath>
#include <iostream>
#include <utils/tiny_tools.h>
#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/dm_mask.h>
#include <Rd/multi_refinable.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  // helper routine to check Alpha
  template <int d, int dt>
  double DKUBasis<d, dt>::alpha_(const int m, const int r) const
  {
    double result(0);

    if (r == 0)
      result = 1.0; // (5.1.1)
    else
      {
	if (m == 0)
	  {
	    double dummy(0);
	    for (int k(ell1()); k <= ell2(); k++)
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
  template <int d, int dt>
  double DKUBasis<d, dt>::alphaT_(const int m, const int r) const
  {
    double result(0);

    if (r == 0)
      result = 1.0; // (5.1.1)
    else
      {
	if (m == 0)
	  {
	    double dummy(0);
	    for (int k(ell1T()); k <= ell2T(); k++)
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
  template <int d, int dt>
  double DKUBasis<d, dt>::betaL_(const int m, const int r) const
  {
    double result(0);

    for (int q((int)ceil((m-ell2T())/2.0)); q <= ellT()-1; q++)
      result += Alpha_(q, r) * cdf_.aT().get_coefficient(m-2*q); // (3.2.31)

    return result * M_SQRT1_2;
  }

  // helper routine to check BetaL
  template <int d, int dt>
  double DKUBasis<d, dt>::betaLT_(const int m, const int r) const
  {
    double result(0);

    for (int q((int)ceil((m-ell2())/2.0)); q <= ell()-1; q++)
      result += AlphaT_(q, r) * cdf_.a().get_coefficient(m-2*q); // (3.2.31)

    return result * M_SQRT1_2;
  }

  template <int d, int dt>
  DKUBasis<d, dt>::DKUBasis(const Array1D<int>& primalBC,
			    const Array1D<int>& dualBC)
    : primalBC_(primalBC),
      dualBC_(dualBC)
  {
    ell1_ = -(d/2);
    ell2_ = d-d/2;
    ell1T_ = ell1_-dt+1;
    ell2T_ = ell2_+dt-1;
    ellT_ = ell2T_;
    ell_ = ellT_-(dt-d);
    
    // setup Alpha
    Alpha_.resize(ellT()+ell2T()-1, dt);
    for (unsigned int m(0); m < Alpha_.row_dimension(); m++)
      Alpha_(m, 0) = 1.0; // (5.1.1)
    for (int r(1); r < (int)Alpha_.column_dimension(); r++)
      {
	double dummy(0);
	for (int k(ell1()); k <= ell2(); k++)
	  {
	    double dummy1(0);
	    for (int s(0); s <= r-1; s++)
	      dummy1 += binomial(r, s) * intpower(k, r-s) * Alpha_(ell2T()-1, s);
	    dummy += cdf_.a().get_coefficient(k) * dummy1; // (5.1.3)
	  }
	Alpha_(ell2T()-1, r) = dummy / (ldexp(1.0, r+1) - 2.0);
      }
    for (int r(1); r < (int)Alpha_.column_dimension(); r++)
      {
	for (int m(-ell2T()+1); m < 0; m++)
	  {
	    double dummy(0);
	    for (int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * Alpha_(ell2T()-1, r-i); // (5.1.2)
	    Alpha_(m+ell2T()-1, r) = dummy;
	  }
	for (int m(1); m <= ellT()-1; m++)
	  {
	    double dummy(0);
	    for (int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * Alpha_(ell2T()-1, r-i); // (5.1.2)
	    Alpha_(m+ell2T()-1, r) = dummy;
	  }
      }

    // setup AlphaT
    AlphaT_.resize(ell()+ell2()-1, d);
    for (unsigned int m(0); m < AlphaT_.row_dimension(); m++)
      AlphaT_(m, 0) = 1.0; // (5.1.1)
    for (int r(1); r < (int)AlphaT_.column_dimension(); r++)
      {
	double dummy(0);
	for (int k(ell1T()); k <= ell2T(); k++)
	  {
	    double dummy1(0);
	    for (int s(0); s <= r-1; s++)
	      dummy1 += binomial(r, s) * intpower(k, r-s) * AlphaT_(ell2()-1, s); // (5.1.3)
	    dummy += cdf_.aT().get_coefficient(k) * dummy1;
	  }
	AlphaT_(ell2()-1, r) = dummy / (ldexp(1.0, r+1) - 2.0);
      }
    for (int r(1); r < (int)AlphaT_.column_dimension(); r++)
      {
	for (int m(-ell2()+1); m < 0; m++)
	  {
	    double dummy(0);
	    for (int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * AlphaT_(ell2()-1, r-i); // (5.1.2)
	    AlphaT_(m+ell2()-1, r) = dummy;
	  }
	for (int m(1); m <= ell()-1; m++)
	  {
	    double dummy(0);
	    for (int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * AlphaT_(ell2()-1, r-i); // (5.1.2)
	    AlphaT_(m+ell2()-1, r) = dummy;
	  }
      }
    
    // setup BetaL
    BetaL_.resize(ell2T()-ell1T()-1, dt);
    for (int r(0); r < dt; r++)
      for (int m(2*ellT()+ell1T()); m <= 2*ellT()+ell2T()-2; m++)
	BetaL_(m-2*ellT()-ell1T(), r) = betaL_(m, r);

    // setup BetaLT
    BetaLT_.resize(ell2()-ell1()-1, d);
    for (int r(0); r < d; r++)
      for (int m(2*ell()+ell1()); m <= 2*ell()+ell2()-2; m++)
	BetaLT_(m-2*ell()-ell1(), r) = betaLT_(m, r);
    
    // setup GammaL:
    //
    // 1. compute the integrals
    //      z(s,t) = \int_0^1\phi(x-s)\tilde\phi(x-t)\,dx
    //    exactly with the [DM] trick
    MultivariateRefinableFunction<DMMask2<HaarMask, CDFMask_primal<d>, CDFMask_dual<d, dt> >, 2> zmask;
    InfiniteVector<double, MultiIndex<int, 2> > zvalues(zmask.evaluate());
    //
    // 2. compute the integrals
    //      I(nu,mu) = \int_0^\infty\phi(x-\nu)\tilde\phi(x-\mu)\,dx
    //    exactly using the z(s,t) values
    Matrix<double> I(ell()-d+dt-1+ell2(), ellT()-1+ell2T());
    for (int nu(-ell2()+1); nu <= -ell1()-1; nu++)
      for (int mu(-ell2T()+1); mu <= -ell1T()-1; mu++)
 	{
 	  double help(0);
 	  int diff(mu - nu);
 	  for (int s(max(-ell2()+1, -ell2T()+1-diff)); s <= nu; s++)
 	    help += zvalues.get_coefficient(MultiIndex<int, 2>(s, s+diff));
	  
 	  I(nu+ell2()-1, mu+ell2T()-1) = help; // (5.1.7)
 	}
    for (int nu(-ell2()+1); nu <= ell()-d+dt-1; nu++)
      for (int mu(-ell2T()+1); mu <= ellT()-1; mu++)
	{
	  if ((nu >= -ell1()) || ((nu <= ell()-1) && (mu >= -ell1T())))
	    I(nu+ell2()-1, mu+ell2T()-1) = (nu == mu ? 1 : 0);
	}
    //
    // 3. 
    
    
//     // check Alpha
//     cout << Alpha_ << endl;
//     for (int m(-ell2T()+1); m <= ellT()-1; m++) {
//       for (int r(0); r < dt; r++) {
// 	cout << alpha_(m, r) << " ";
//       }
//       cout << endl;
//     }
    
//     cout << AlphaT_ << endl;
//     for (int m(-ell2()+1); m <= ell()-1; m++) {
//       for (int r(0); r < d; r++) {
// 	cout << alphaT_(m, r) << " ";
//       }
//       cout << endl;
//     }

//     // check BetaL
//     cout << BetaL_ << endl;

//     // check BetaLT
//     cout << BetaLT_ << endl;
  }
}

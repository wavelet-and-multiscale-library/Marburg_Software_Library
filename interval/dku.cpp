// implementation for dku.h

#include <cmath>
#include <iostream>
#include <utils/tiny_tools.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  // helper routine to check Alpha
  template <unsigned int d, unsigned int dt>
  double DKUBasis<d, dt>::alpha_(const int m, const unsigned int r) const
  {
    double result(0);

    if (r == 0u)
      result = 1.0; // (5.1.1)
    else
      {
	if (m == 0)
	  {
	    double dummy(0);
	    for (int k(ell1()); k <= ell2(); k++)
	      {
		double dummy1(0);
		for (unsigned int s(0); s <= r-1; s++)
		  dummy1 += binomial(r, s) * intpower(k, r-s) * alpha_(0, s);
		dummy += cdf_.a().get_coefficient(k) * dummy1;
	      }
	    result = dummy / (ldexp(1.0, r+1) - 2.0);
	  }
	else
	  {
 	    for (unsigned int i(0); i <= r; i++)
	      result += binomial(r, i) * intpower(m, i) * alpha_(0, r-i); // (5.1.2)
	  }
      }

    return result;
  }

   // helper routine to check AlphaT
  template <unsigned int d, unsigned int dt>
  double DKUBasis<d, dt>::alphaT_(const int m, const unsigned int r) const
  {
    double result(0);

    if (r == 0u)
      result = 1.0; // (5.1.1)
    else
      {
	if (m == 0)
	  {
	    double dummy(0);
	    for (int k(ell1T()); k <= ell2T(); k++)
	      {
		double dummy1(0);
		for (unsigned int s(0); s <= r-1; s++)
		  dummy1 += binomial(r, s) * intpower(k, r-s) * alphaT_(0, s);
		dummy += cdf_.aT().get_coefficient(k) * dummy1;
	      }
	    result = dummy / (ldexp(1.0, r+1) - 2.0);
	  }
	else
	  {
 	    for (unsigned int i(0); i <= r; i++)
	      result += binomial(r, i) * intpower(m, i) * alphaT_(0, r-i); // (5.1.2)
	  }
      }

    return result;
  }

 template <unsigned int d, unsigned int dt>
  DKUBasis<d, dt>::DKUBasis()
    : ellT_(ell2T()),
      Alpha_(ellT()+ell2T()-1, dt),
      AlphaT_(ell()+ell2()-1, d)
  {
    // setup Alpha
    for (unsigned int m(0); m < Alpha_.row_dimension(); m++)
      Alpha_(m, 0) = 1.0; // (5.1.1)
    for (unsigned int r(1); r < Alpha_.column_dimension(); r++)
      {
	double dummy(0);
	for (int k(ell1()); k <= ell2(); k++)
	  {
	    double dummy1(0);
	    for (unsigned int s(0); s <= r-1; s++)
	      dummy1 += binomial(r, s) * intpower(k, r-s) * Alpha_(ell2T()-1, s); // (5.1.3)
	    dummy += cdf_.a().get_coefficient(k) * dummy1;
	  }
	Alpha_(ell2T()-1, r) = dummy / (ldexp(1.0, r+1) - 2.0);
      }
    for (unsigned int r(1); r < Alpha_.column_dimension(); r++)
      {
	for (int m(-ell2T()+1); m < 0; m++)
	  {
	    double dummy(0);
	    for (unsigned int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * Alpha_(ell2T()-1, r-i); // (5.1.2)
	    Alpha_(m+ell2T()-1, r) = dummy;
	  }
	for (int m(1); m <= ellT()-1; m++)
	  {
	    double dummy(0);
	    for (unsigned int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * Alpha_(ell2T()-1, r-i); // (5.1.2)
	    Alpha_(m+ell2T()-1, r) = dummy;
	  }
      }

    // setup AlphaT
    for (unsigned int m(0); m < AlphaT_.row_dimension(); m++)
      AlphaT_(m, 0) = 1.0; // (5.1.1)
    for (unsigned int r(1); r < AlphaT_.column_dimension(); r++)
      {
	double dummy(0);
	for (int k(ell1T()); k <= ell2T(); k++)
	  {
	    double dummy1(0);
	    for (unsigned int s(0); s <= r-1; s++)
	      dummy1 += binomial(r, s) * intpower(k, r-s) * AlphaT_(ell2()-1, s); // (5.1.3)
	    dummy += cdf_.aT().get_coefficient(k) * dummy1;
	  }
	AlphaT_(ell2()-1, r) = dummy / (ldexp(1.0, r+1) - 2.0);
      }
    for (unsigned int r(1); r < AlphaT_.column_dimension(); r++)
      {
	for (int m(-ell2()+1); m < 0; m++)
	  {
	    double dummy(0);
	    for (unsigned int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * AlphaT_(ell2()-1, r-i); // (5.1.2)
	    AlphaT_(m+ell2()-1, r) = dummy;
	  }
	for (int m(1); m <= ell()-1; m++)
	  {
	    double dummy(0);
	    for (unsigned int i(0); i <= r; i++)
	      dummy += binomial(r, i) * intpower(m, i) * AlphaT_(ell2()-1, r-i); // (5.1.2)
	    AlphaT_(m+ell2()-1, r) = dummy;
	  }
      }

    // check Alpha
    cout << Alpha_ << endl;
    for (int m(-ell2T()+1); m <= ellT()-1; m++) {
      for (unsigned int r(0); r < dt; r++) {
	cout << alpha_(m, r) << " ";
      }
      cout << endl;
    }
    
    cout << AlphaT_ << endl;
    for (int m(-ell2()+1); m <= ell()-1; m++) {
      for (unsigned int r(0); r < d; r++) {
	cout << alphaT_(m, r) << " ";
      }
      cout << endl;
    }
  }
}

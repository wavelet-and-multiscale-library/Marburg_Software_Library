// implementation for dku.h

#include <cmath>
#include <iostream>
#include <utils/tiny_tools.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  // just a helper routine to check Alpha
  template <unsigned int d, unsigned int dt>
  double DKUBasis<d, dt>::Alpha_(const int m, const unsigned int r)
  {
    if (r == 0u)
      return 0; // (5.1.1)
    else
      {
	
      }

    return 42;
  }

  template <unsigned int d, unsigned int dt>
  DKUBasis<d, dt>::DKUBasis()
    : ellT_(ell2T()),
      Alpha_(ellT()+ell2T()-1, dt),
      AlphaT_(ell()+ell2()-1, d)
  {
    // setup the moment matrix Alpha
    for (unsigned int m(0); m < Alpha_.row_dimension(); m++)
      Alpha_(m, 0) = 1.0; // (5.1.1)
    for (unsigned int r(1); r < Alpha_.column_dimension(); r++)
      {
	double factor(1.0 / (ldexp(1.0, r+1) - 2.0));
	double dummy(0);
	for (int k(ell1()); k <= ell2(); k++)
	  {
	    double dummy1(0);
	    for (unsigned int s(0); s <= r-1; s++)
	      dummy1 += binomial(r, s) * intpower(k, (int)r-(int)s) * Alpha_(ell2T()-1, s); // (5.1.3)
	    dummy += cdf_.a().get_coefficient(k) * dummy1;
	  }
	Alpha_(ell2T()-1, r) = factor * dummy;
      }

    cout << Alpha_ << endl;

    for (int m(-ell2T()+1); m <= ellT()-1; m++)
      {
	for (unsigned int r(0); r < dt; r++)
	  {
	    cout << Alpha_(m, r) << " ";
	  }
	cout << endl;
      }
    
  }
}

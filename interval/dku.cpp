// implementation for dku.h

#include <cmath>
#include <utils/tiny_tools.h>

namespace WaveletTL
{
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
  }
}

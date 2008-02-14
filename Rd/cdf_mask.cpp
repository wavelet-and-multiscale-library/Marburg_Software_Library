// implementation for cdf_mask.h

#include <cassert>
#include <cmath>
#include <Rd/cdf_utils.h>
#include <algebra/polynomial.h>
#include <utils/tiny_tools.h>
#include <utils/multiindex.h>

using MathTL::Polynomial;
using MathTL::MultiIndex;

namespace WaveletTL
{
  template <int d>
  CDFMask_primal<d>::CDFMask_primal()
  {
    // The primal generators are B-splines of order d, centered around (d mod 2)/2.
    // So in z-notation, the mask of the primal generators is just
    //   a(z) = 2^{1-d} * (1+z)^d * z^{-floor(d/2)},
    // the coefficients of which can be calculated in closed form (faster):
    for (int k(ell1<d>()); k <= ell2<d>(); k++)
      set_coefficient(MultiIndex<int, 1>(k), ldexp(1., 1 - d) * binomial(d, k-ell1<d>()));
  }
  
  template <int d>
  CDFRefinementMask_primal<d>::CDFRefinementMask_primal()
  {
    // The primal generators are B-splines of order d, centered around (d mod 2)/2.
    // So in z-notation, the mask of the primal generators is just
    //   a(z) = (1+z)^d * z^{-floor(d/2)},
    // the coefficients of which can be calculated in closed form (faster):
    for (int k(ell1<d>()); k <= ell2<d>(); k++)
      RRefinementMask<d+1,-(d/2)>::coeffs_[k-ell1<d>()] = ldexp(1., 1 - d) * binomial(d, k-ell1<d>());
  }
  
  template <int d, int dt>
  CDFMask_dual<d, dt>::CDFMask_dual()
  {
    assert(d >= 1);
    assert(dt >= 0);
    assert((d%2) == (dt%2));
    
    // Construct the mask of the dual generator, cf. [CDF] for details.
#if 0
    // (old version)
    Polynomial<double> u, uh1, uh2, uh2a(1), uh3;
    
    // uh1 = z
    uh1.set_coefficient(1, 1);
    
    // uh2 = (z-1)^2
    uh2.set_coefficient(0, -1);
    uh2.set_coefficient(1, 1);
    uh2.multiply(uh2);
    
    // uh3 = z+1
    uh3.set_coefficient(0, 1);
    uh3.set_coefficient(1, 1);
    
    int K = (d+dt)/2;
    
    for (int k(0); k<K; k++)
      {
 	u.add((binomial(K-1+k, k) / pow(-4.0, k)), uh1.power(K-1-k) * uh2a);
 	uh2a.multiply(uh2);
      }
    
    u.multiply(ldexp(1.0, -dt+1) * uh3.power(dt));
    
    for (unsigned int k(0); k <= u.degree(); k++)
      set_coefficient(MultiIndex<int, 1>((int)k - d/2 - dt + 1), u.get_coefficient(k));
#else
    // use explicit formula from [P]
    for (int k(ell1T<d,dt>()); k <= ell2T<d,dt>(); k++) {
      double help = 0;
      for (int n = 0; n <= (d+dt)/2-1; n++)
	for (int ell = 0; ell <= 2*n; ell++)
	  help +=
	    ldexp(1.0, 1-dt-2*n) // typo in [P]
	    * minus1power(n+ell)
	    * binomial(dt, k+dt/2-ell+n)
	    * binomial((d+dt)/2-1+n, n)
	    * binomial(2*n, ell);
      
      set_coefficient(MultiIndex<int, 1>(k), help);
    }
#endif
  }

  template <int d, int dt>
  CDFRefinementMask_dual<d, dt>::CDFRefinementMask_dual()
  {
    assert(d >= 1);
    assert(dt >= 0);
    assert((d%2) == (dt%2));
    
    // Construct the mask of the dual generator, cf. [CDF] for details.
    // Here we use the explicit formula from [P]
    for (int k(ell1T<d,dt>()); k <= ell2T<d,dt>(); k++) {
      double help = 0;
      for (int n = 0; n <= (d+dt)/2-1; n++)
	for (int ell = 0; ell <= 2*n; ell++)
	  help +=
	    ldexp(1.0, 1-dt-2*n) // typo in [P]
	    * minus1power(n+ell)
	    * binomial(dt, k+dt/2-ell+n)
	    * binomial((d+dt)/2-1+n, n)
	    * binomial(2*n, ell);
      
      RRefinementMask<d+2*dt-1,-(d/2)-dt+1>::coeffs_[k-ell1T<d,dt>()] = help;
    }
  }
  
}

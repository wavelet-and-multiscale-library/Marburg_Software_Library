// implementation for spline_expansion.h

#include <numerics/quadrature.h>
#include <numerics/iteratsolv.h>
#include <galerkin/full_gramian.h>

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,flavor>& basis,
	 const bool primal,
	 const int jmax,
	 Vector<double>& coeffs)
  {
    basis.expand(f, primal, jmax, coeffs);
  }

  template <int d, int dT, SplineBasisFlavor flavor>
  inline
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,flavor>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs)
  {
    basis.expand(f, primal, jmax, coeffs);
  }
  
}

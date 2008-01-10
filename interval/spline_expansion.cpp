// implementation for spline_expansion.h

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
	 const bool primal,
	 const int jmax,
	 Vector<double>& coeffs)
  {
    basis.expand(f, primal, jmax, coeffs);
  }

  template <int d, int dT, SplineBasisFlavor flavor, int s0, int s1, int sT0, int sT1>
  inline
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename SplineBasis<d,dT,flavor,s0,s1,sT0,sT1>::Index>& coeffs)
  {
    basis.expand(f, primal, jmax, coeffs);
  }
  
}

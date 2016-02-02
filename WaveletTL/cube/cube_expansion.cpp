// implementation for cube_expansion.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  void
  expand(const Function<DIM>* f,
	 const CubeBasis<IBASIS,DIM>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs)
  {
    basis.expand(f, primal, jmax, coeffs);
  }

}

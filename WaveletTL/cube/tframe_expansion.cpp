// implementation for cube_expansion.h

namespace WaveletTL
{
  template <class IFRAME, unsigned int DIM>
  void
  expand(const Function<DIM>* f,
	 const TensorFrame<IFRAME,DIM>& frame,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename TensorFrame<IFRAME,DIM>::Index>& coeffs)
  {
    frame.expand(f, primal, jmax, coeffs);
  }

}

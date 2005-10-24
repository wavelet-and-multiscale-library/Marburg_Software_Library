// implementation for cube_evaluate.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const CubeBasis<IBASIS,DIM>& basis,
			       const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			       const bool primal,
			       const int resolution)
  {
    
    return SampledMapping<DIM>(); // dummy return for the compiler
  }

}

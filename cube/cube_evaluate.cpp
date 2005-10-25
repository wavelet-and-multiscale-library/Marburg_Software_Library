// implementation for cube_evaluate.h

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const CubeBasis<IBASIS,DIM>& basis,
			       const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			       const bool primal,
			       const int resolution)
  {
    Grid<DIM> grid(Point<DIM>(0), Point<DIM>(1), 1<<resolution);
//     Grid<1> grid(0, 1, 1<<resolution);
//     Array1D<double> values((1<<resolution)+1);

  
    return SampledMapping<DIM>(); // dummy return for the compiler
  }

}

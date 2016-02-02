// implementation for cube_evaluate.h

#include <utils/array1d.h>
#include <utils/fixed_array1d.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

using namespace MathTL;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda
    for (unsigned int i = 0; i < DIM; i++)
      values[i] = evaluate(*(basis.bases()[i]),
			   typename IBASIS::Index(lambda.j(),
						  lambda.e()[i],
						  lambda.k()[i],
						  basis.bases()[i]),
			   primal,
			   resolution).values();
      
    SampledMapping<DIM> result(Point<DIM>(0), Point<DIM>(1), values);
    return result; // gcc 2.95 does not like these two lines melted into one
  }
    
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM>
  evaluate(const CubeBasis<IBASIS,DIM>& basis,
	   const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    Grid<DIM> grid(Point<DIM>(0), Point<DIM>(1), 1<<resolution);
    SampledMapping<DIM> result(grid); // zero
      
    typedef typename CubeBasis<IBASIS,DIM>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	   itend(coeffs.end()); it != itend; ++it)
      result.add(*it, evaluate(basis, it.index(), primal, resolution));
      
    return result;
  }
}

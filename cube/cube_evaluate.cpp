// implementation for cube_evaluate.h

#include <utils/array1d.h>
#include <utils/fixed_array1d.h>

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const CubeBasis<IBASIS,DIM>& basis,
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
    
    return SampledMapping<DIM>(Point<DIM>(0), Point<DIM>(1), values);
  }

}

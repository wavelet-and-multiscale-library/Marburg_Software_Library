// implementation for frame_evaluate.h

#include <utils/array1d.h>
#include <utils/fixed_array1d.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

namespace FrameTL
{

 template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
 SampledMapping<DIM_d>
 evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
	  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
	  const bool primal,
	  const int resolution)
 {
   assert(DIM_d == 2  || DIM_d == 1);
   
   FixedArray1D<Array1D<double>,DIM_d> values; // point values of the factors within psi_lambda
   for (unsigned int i = 0; i < DIM_d; i++)
     values[i] = evaluate(*(frame.bases()[lambda.p()]->bases()[i]),
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame.bases()[lambda.p()]->bases()[i]),
			  primal,
			  resolution).values();   

   SampledMapping<DIM_d> result(*(frame.atlas()->charts()[lambda.p()]),
				values,
				resolution);
   
   return result; // gcc 2.95 does not like these two lines melted into one
 }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  Array1D<SampledMapping<DIM_d> >
  evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
	   const InfiniteVector<double,
	   typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {

    Array1D<SampledMapping<DIM_d> > result(frame.n_p());

    for (unsigned int i = 0; i < frame.n_p(); i++) {
      result[i] = SampledMapping<DIM_d>(*(frame.atlas()->charts()[i]),resolution);// all zero
    }

    typedef typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
 	  itend(coeffs.end()); it != itend; ++it) {
      result[it.index().p()].add(*it, evaluate(frame, it.index(), primal, resolution));
    }

    return result;
  }

}

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
	  const unsigned int patch,
	  const bool primal,
	  const int resolution)
 {
   assert(DIM_d == 2  || DIM_d == 1);
   
   //SampledMapping<DIM_d> res(*(frame.atlas()->charts()[patch]),resolution);// all zero
   switch (DIM_d) {
   case 2: {
     Matrix<double> gridx;
     Matrix<double> gridy;
    // setup grid
     const unsigned int  n_points = (1 << resolution)+1;
     const double h = 1.0 / (n_points-1);
     
     Point<DIM_d> x;
     Point<DIM_m> x_patch;
     Point<DIM_d> y;

     Matrix<double> values;
     values.resize(n_points,n_points);

     gridx.resize(n_points,n_points);
     gridy.resize(n_points,n_points);

     for (unsigned int i = 0; i < n_points; i++) {
       x[0] = h*i;
       for (unsigned int j = 0; j < n_points; j++) {
	 x[1] = h*j;
	 frame.atlas()->charts()[patch]->map_point(x,x_patch);
	 gridx.set_entry(i,j,x_patch[0]);
	 gridy.set_entry(i,j,x_patch[1]);
	 if ( in_support(frame, lambda, x_patch) ) {
	   frame.atlas()->charts()[lambda.p()]->map_point_inv(x_patch,y);
	   double wav_val_x = WaveletTL::evaluate(*(frame.bases()[lambda.p()]->bases()[0]), 0,
						  typename IBASIS::Index(lambda.j(),
									 lambda.e()[0],
									 lambda.k()[0],
									 frame.bases()[lambda.p()]->bases()[0]),
						  y[0]);
	   double wav_val_y = WaveletTL::evaluate(*(frame.bases()[lambda.p()]->bases()[1]), 0,
						  typename IBASIS::Index(lambda.j(),
									 lambda.e()[1],
									 lambda.k()[1],
									 frame.bases()[lambda.p()]->bases()[1]),
						  y[1]);
	   values.set_entry(i,j,(wav_val_x*wav_val_y) / frame.atlas()->charts()[lambda.p()]->Gram_factor(y));
	 }
	 else {
	   values.set_entry(i,j,0.);
	 }
       }
     }
     return SampledMapping<2> (Grid<2>(gridx,gridy),values);
   }// end case 2
     
   case 1: {
     //TODO
   }
   }


 }

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  SampledMapping<DIM_d> evaluate_single_patch(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				 const InfiniteVector<double,
				 typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& coeffs,
				 const bool primal,
				 const int resolution)
  {
    assert(DIM_d == 2  || DIM_d == 1);
    typedef typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index Index;
    typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
      itend(coeffs.end());
    
    SampledMapping<DIM_d> result(*(frame.atlas()->charts()[it.index().p()]),resolution);// all zero
    
    for (; it != itend; ++it)
      result.add(*it, evaluate(frame, it.index(), primal, resolution));

    return result;
  }


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
      for (unsigned int i = 0; i < frame.n_p(); i++) {
	if ( frame.atlas()->get_adjacency_matrix().get_entry(i,it.index().p()) ) {
	  if ( i == it.index().p())
	    result[i].add(*it, evaluate(frame, it.index(), primal, resolution));
	  else
	    result[i].add(*it, evaluate(frame, it.index(), i, primal, resolution));
	}
      }
    }

    return result;
  }

}

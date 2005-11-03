// implementation for frame_evaluate.h

#include <utils/array1d.h>
#include <utils/fixed_array1d.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

namespace FrameTL
{

  template <class IBASIS>
  SampledMapping<1>
  EvaluateFrame<IBASIS,1,1> ::evaluate(const AggregatedFrame<IBASIS,1,1>& frame,
				       const typename AggregatedFrame<IBASIS,1,1>::Index& lambda,
				       const unsigned int patch,
				       const bool primal,
				       const int resolution) const
  {
    
    Point<1> x;
    Point<1> x_patch;
    Point<1> y;

    const unsigned int  n_points = (1 << resolution)+1;
    const double h = 1.0 / (n_points-1);
    Array1D<double> grid(n_points);
    Array1D<double> values(n_points);
    for (unsigned int i = 0; i < n_points; i++) {
      x[0] = h*i;
      frame.atlas()->charts()[patch]->map_point(x,x_patch);
      if ( in_support(frame, lambda, x_patch) ) {
	frame.atlas()->charts()[lambda.p()]->map_point_inv(x_patch,y);
	double wav_val = WaveletTL::evaluate(*(frame.bases()[lambda.p()]->bases()[0]), 0,
					     typename IBASIS::Index(lambda.j(),
								    lambda.e()[0],
								    lambda.k()[0],
								    frame.bases()[lambda.p()]->bases()[0]),
					     
					     y[0]);
	values[i] = wav_val / frame.atlas()->charts()[lambda.p()]->Gram_factor(y);
      }
      else
	values[i] = 0.;
    }
    
    return SampledMapping<1> (Grid<1>(grid),values);
  }
  

 template <class IBASIS>
 SampledMapping<1>
 EvaluateFrame<IBASIS,1,1>::evaluate(const AggregatedFrame<IBASIS,1,1>& frame,
					     const typename AggregatedFrame<IBASIS,1,1>::Index& lambda,
					     const bool primal,
					     const int resolution) const
 {
 
   FixedArray1D<Array1D<double>,1> values; // point values of the factors within psi_lambda
   for (unsigned int i = 0; i < 1; i++)
     values[i] = WaveletTL::evaluate(*(frame.bases()[lambda.p()]->bases()[i]),
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame.bases()[lambda.p()]->bases()[i]),
			  primal,
			  resolution).values();

   SampledMapping<1> result(*(frame.atlas()->charts()[lambda.p()]),
				values,
				resolution);
   
   return result; // gcc 2.95 does not like these two lines melted into one
 }


  template <class IBASIS>
  Array1D<SampledMapping<1> >
  EvaluateFrame<IBASIS,1,1>::evaluate(const AggregatedFrame<IBASIS,1,1>& frame,
	   const InfiniteVector<double,
	   typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
	   const bool primal,
	   const int resolution) const
  {

    Array1D<SampledMapping<1> > result(frame.n_p());

    for (unsigned int i = 0; i < frame.n_p(); i++) {      
      result[i] = SampledMapping<1>(*(frame.atlas()->charts()[i]),resolution);// all zero
    }

    typedef typename AggregatedFrame<IBASIS,1,1>::Index Index;
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

  
  template <class IBASIS>
  SampledMapping<2>
  EvaluateFrame<IBASIS,2,2> ::evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
				       const typename AggregatedFrame<IBASIS,2,2>::Index& lambda,
				       const unsigned int patch,
				       const bool primal,
				       const int resolution) const
  {

    Point<2> x;
    Point<2> x_patch;
    Point<2> y;

    const unsigned int  n_points = (1 << resolution)+1;
    const double h = 1.0 / (n_points-1);     
    
    Matrix<double> gridx;
    Matrix<double> gridy;
    // setup grid
    
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
  }

 template <class IBASIS>
 SampledMapping<2>
 EvaluateFrame<IBASIS,2,2>::evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
				     const typename AggregatedFrame<IBASIS,2,2>::Index& lambda,
				     const bool primal,
				     const int resolution) const
 {
   FixedArray1D<Array1D<double>,2> values; // point values of the factors within psi_lambda
   for (unsigned int i = 0; i < 2; i++)
     values[i] = WaveletTL::evaluate(*(frame.bases()[lambda.p()]->bases()[i]),
			  typename IBASIS::Index(lambda.j(),
						 lambda.e()[i],
						 lambda.k()[i],
						 frame.bases()[lambda.p()]->bases()[i]),
			  primal,
			  resolution).values();

   SampledMapping<2> result(*(frame.atlas()->charts()[lambda.p()]),
				values,
				resolution);
   
   return result; // gcc 2.95 does not like these two lines melted into one
 }

  template <class IBASIS>
  Array1D<SampledMapping<2> >
  EvaluateFrame<IBASIS,2,2>::evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
	   const InfiniteVector<double,
	   typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
	   const bool primal,
	   const int resolution) const
  {

    Array1D<SampledMapping<2> > result(frame.n_p());

    for (unsigned int i = 0; i < frame.n_p(); i++) {      
      result[i] = SampledMapping<2>(*(frame.atlas()->charts()[i]),resolution);// all zero
    }

    typedef typename AggregatedFrame<IBASIS,2,2>::Index Index;
    for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
 	  itend(coeffs.end()); it != itend; ++it) {
     
      for (unsigned int i = 0; i < frame.n_p(); i++) {
	if ( frame.atlas()->get_adjacency_matrix().get_entry(i,it.index().p()) ) {
	  if ( i == it.index().p()) {
	    result[i].add(*it, evaluate(frame, it.index(), primal, resolution));
	  }
	  else
	    result[i].add(*it, evaluate(frame, it.index(), i, primal, resolution));
	}
      }
    }

    return result;
  }

}

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_FRAME_EVALUATE_H
#define _FRAMETL_FRAME_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <aggregated_frame.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;
using FrameTL::AggregatedFrame;

namespace FrameTL
{

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class EvaluateFrame
  {

  public:    
    /*!
      Evaluate a single primal/dual generator or wavelet \psi_lambda
      on a dyadic grid of the patch given by 'patch'.
    */
    virtual SampledMapping<DIM_d> evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
					   const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
					   const unsigned int patch,
					   const bool primal,
					   const int resolution) const = 0;
    
    /*!
      Evaluate a single primal/dual generator or wavelet \psi_\lambda
      on a dyadic subgrid of its corresponding patch.
    */
    virtual SampledMapping<DIM_d> evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
					   const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
					   const bool primal,
					   const int resolution) const = 0;
    
    /*!
      Evaluate an arbitrary linear combination of primal/dual wavelets
      on a dyadic subgrid of each patch.
    */
    virtual Array1D<SampledMapping<DIM_d> > evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
						     const InfiniteVector<double,
						     typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& coeffs,
						     const bool primal,
						     const int resolution) const = 0;

    /*!
      Evaluates the difference between the function given by the expansion coefficients
      'coeffs' and the function f
    */
    virtual Array1D<SampledMapping<DIM_d> > evaluate_difference(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
								const InfiniteVector<double,
								typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& coeffs,
								const Function<DIM_m>& f,
								const int resolution) const = 0;


    virtual double  L_2_error(const AggregatedFrame<IBASIS,1,1>& frame,
			      const InfiniteVector<double,
			      typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
			      const Function<1>& f,
			      const int resolution,
			      const double a, const double b) const;
  };

  template <class IBASIS>
  class EvaluateFrame<IBASIS,1,1>
  {
  public:

    SampledMapping<1>
    evaluate(const AggregatedFrame<IBASIS,1,1>& frame,
	     const typename AggregatedFrame<IBASIS,1,1>::Index& lambda,
	     const unsigned int patch,
	     const bool primal,
	     const int resolution) const;

    SampledMapping<1>
    evaluate(const AggregatedFrame<IBASIS,1,1>& frame,
	     const typename AggregatedFrame<IBASIS,1,1>::Index& lambda,
	     const bool primal,
	     const int resolution) const;
    
    Array1D<SampledMapping<1> >
    evaluate(const AggregatedFrame<IBASIS,1,1>& frame,
	     const InfiniteVector<double,
	     typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
	     const bool primal,
	     const int resolution) const;

    Array1D<SampledMapping<1> >
    evaluate_difference(const AggregatedFrame<IBASIS,1,1>& frame,
			const InfiniteVector<double,
			typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
			const Function<1>& f,
			const int resolution) const;

    double
    L_2_error(const AggregatedFrame<IBASIS,1,1>& frame,
	      const InfiniteVector<double,
	      typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
	      const Function<1>& f,
	      const int resolution,
	      const double a, const double b) const;



  };
  
  template <class IBASIS>
  class EvaluateFrame<IBASIS,2,2>
  {
  public:

    SampledMapping<2>
    evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
	     const typename AggregatedFrame<IBASIS,2,2>::Index& lambda,
	     const unsigned int patch,
	     const bool primal,
	     const int resolution) const;

    SampledMapping<2>
    evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
	     const typename AggregatedFrame<IBASIS,2,2>::Index& lambda,
	     const bool primal,
	     const int resolution) const;
    
    Array1D<SampledMapping<2> >
    evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
	     const InfiniteVector<double,
	     typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
	     const bool primal,
	     const int resolution) const;

    Array1D<SampledMapping<2> >
    evaluate_difference(const AggregatedFrame<IBASIS,2,2>& frame,
			const InfiniteVector<double,
			typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
			const Function<2>& f,
			const int resolution) const;


    double
    L_2_error(const AggregatedFrame<IBASIS,2,2>& frame,
	      const InfiniteVector<double,
	      typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
	      const Function<2>& f,
	      const int resolution,
	      const double a, const double b) const;

  };

}

#include <frame_evaluate.cpp>

#endif

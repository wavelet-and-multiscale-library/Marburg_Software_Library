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
  /*!
    Abstract base class for classes providing functions for the point
    evaluation of single frame elements or linear combinations of frame elements.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class EvaluateFrame
  {

  public:
    /*!
      Evaluate a single primal/dual generator or wavelet \f$\psi_\lambda\f$
      on a dyadic grid of the patch given by 'patch'.
    */
    virtual SampledMapping<DIM_d> evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
					   const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
					   const unsigned int patch,
					   const bool primal,
					   const int resolution) const = 0;

    /*!
      Evaluate a single primal/dual generator or wavelet \f$\psi_\lambda\f$
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
						     const InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& coeffs,
						     const bool primal,
						     const int resolution) const = 0;

    /*!
      Evaluates the difference between the function given by the expansion coefficients
      'coeffs' and the function f
    */
    virtual Array1D<SampledMapping<DIM_d> > evaluate_difference(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
								const InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& coeffs,
								const Function<DIM_m>& f,
								const int resolution) const = 0;


    virtual double  L_2_error(const AggregatedFrame<IBASIS,1,1>& frame,
			      const InfiniteVector<double, typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
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
	     const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
	     const bool primal,
	     const int resolution) const;

    Array1D<SampledMapping<2> >
    evaluate_difference(const AggregatedFrame<IBASIS,2,2>& frame,
			const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
			const Function<2>& f,
			const int resolution) const;


    /*! C.Hartmann: added on 15.04.16
    Evaluate a single primal frame \psi_\lambda
    on an arbitrary grid. The grid is given by a tensor product
    of 1-dimensional grids. Makes use of tensor product structure of the frames.
    ATTENTION: USE ONLY FOR AGGREGATED FRAME WHICH CHARTS ARE *SIMPLE AFFINE LINEAR*
    */
    Matrix<double>
    evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
			const typename AggregatedFrame<IBASIS,2,2>::Index& lambda,
			const FixedArray1D<Array1D<double>,2>& grid) const;


    /*! C.Hartmann: added on 15.04.16
    Evaluate an arbitrary linear combination of primal frames
    on an arbitrary grid. The grid is given by a tensor product
    of 1-dimensional grids. Makes use of tensor product structure of the frames.
    ATTENTION: USE ONLY FOR AGGREGATED FRAME WHICH CHARTS ARE *SIMPLE AFFINE LINEAR*
    */
    Matrix<double>
    evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
             const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
             const FixedArray1D<Array1D<double>,2>& grid) const;


    /*! C.Hartmann: added on 15.04.16
    Evaluate the first derivative of an arbitrary linear combination of primal frames
    on an arbitrary grid. The grid is given by a tensor product
    of 1-dimensional grids. Makes use of tensor product structure of the frames.
    ATTENTION: USE ONLY FOR AGGREGATED FRAME WHICH CHARTS ARE *SIMPLE AFFINE LINEAR*
    */
    void
    evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
             const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
             const FixedArray1D<Array1D<double>,2>& grid,
             Matrix<double>& deriv_x, Matrix<double>& deriv_y) const;


    double
    L_2_error(const AggregatedFrame<IBASIS,2,2>& frame,
	      const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
	      const Function<2>& f,
	      const int resolution,
	      const double a, const double b) const;


    double
    evaluate(const AggregatedFrame<IBASIS,2,2>& frame,
			const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
			const Point<2> x) const;


    void
    evaluate_deriv(const AggregatedFrame<IBASIS,2,2>& frame,
			const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
			const Point<2> x,
			FixedArray1D<double,2>& deriv_values) const;

    void
    evaluate_deriv_tuned(const AggregatedFrame<IBASIS,2,2>& frame,
			const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
			const Point<2> x,
			FixedArray1D<double,2>& deriv_values) const;

  };

}

#include <frame_evaluate.cpp>

#endif

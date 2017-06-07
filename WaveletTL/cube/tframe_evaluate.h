// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TFRAME_EVALUATE_H
#define	_WAVELETTL_TFRAME_EVALUATE_H


#include <algebra/infinite_vector.h>
#include <cube/tframe.h>
#include <geometry/sampled_mapping.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <utils/array1d.h>
#include <utils/fixed_array1d.h>
#include <geometry/grid.h>

using MathTL::Grid;
using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <class IFRAME, unsigned int DIM> class TensorFrame;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1]^d.
  */
  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorFrame<IFRAME,DIM>& frame,
			       const typename TensorFrame<IFRAME,DIM>::Index& lambda,
			       const bool primal,
			       const int resolution);
  
  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorFrame<IFRAME,DIM>& frame,
			       const int& lambdanum,
			       const bool primal,
			       const int resolution);

  /*!
    Evaluate an arbitrary linear combination of primal/dual wavelets
    on a dyadic subgrid of [0,1]^d.
  */
  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorFrame<IFRAME,DIM>& frame,
			       const InfiniteVector<double, typename TensorFrame<IFRAME,DIM>::Index>& coeffs,
			       const bool primal,
			       const int resolution);
  
  template <class IFRAME, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorFrame<IFRAME,DIM>& frame,
			       const InfiniteVector<double, int>& coeffs,
			       const bool primal,
			       const int resolution);
}

#include <cube/tframe_evaluate.cpp>

#endif	/* _WAVELETTL_TFRAME_EVALUATE_H */


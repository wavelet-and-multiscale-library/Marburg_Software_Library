// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_EVALUATE_H
#define _WAVELETTL_SPLINE_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT> class SplineBasis;

  /*!
    Evaluate an arbitrary linear combination of primal wavelets
    on a dyadic subgrid of [0,1].
  */
  template <int d, int dT>
  SampledMapping<1> evaluate(const SplineBasis<d,dT>& basis,
			     const InfiniteVector<double, typename SplineBasis<d,dT>::Index>& coeffs,
			     const int resolution);

}

#include <interval/spline_evaluate.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_INTERVAL_EVALUATE_H
#define _WAVELETTL_INTERVAL_EVALUATE_H

#include <utils/array1d.h>
#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;
using MathTL::Array1D;

namespace WaveletTL
{
  /*!
    Evaluate a single primal generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].
  */
  template <class IBasis>
  SampledMapping<1> evaluate(const IBasis& basis,
                             const typename IBasis::Index& lambda,
                             const int resolution);
  
  /*!
    Evaluate a single primal or dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].
  */
  template <class IBasis>
  SampledMapping<1> evaluate(const IBasis& basis,
                             const typename IBasis::Index& lambda,
                             const bool primal,
                             const int resolution);

  /*!
    Evaluate an arbitrary linear combination of primal
    wavelets on a dyadic subgrid of [0,1].
  */
  template <class IBasis>
  SampledMapping<1> evaluate(const IBasis& basis,
                             const InfiniteVector<double, typename IBasis::Index>& coeffs,
                             const int resolution);

  /*!
    dummy routine
  */
  template <class IBasis>
  SampledMapping<1> evaluate(const IBasis& basis,
                             const InfiniteVector<double, typename IBasis::Index>& coeffs,
                             const bool primal,
                             const int resolution)
  {
    return evaluate(basis, coeffs, resolution);
  }

  /*!
    point evaluation of (derivatives) of a single primal generator
    or wavelet \psi_\lambda
  */
  template <class IBasis>
  double evaluate(const IBasis& basis, const unsigned int derivative,
                  const typename IBasis::Index& lambda,
                  const double x);

  /*!
    point evaluation of (derivatives) of a single primal generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <class IBasis>
  void evaluate(const IBasis& basis, const unsigned int derivative,
                const typename IBasis::Index& lambda,
                const Array1D<double>& points, Array1D<double>& values);

  /*!
    point evaluation of 0-th and first derivative of a single primal generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <class IBasis>
  void evaluate(const IBasis& basis,
                const typename IBasis::Index& lambda,
                const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);

}

#include <interval/interval_evaluate.cpp>

#endif // _WAVELETTL_INTERVAL_EVALUATE_H

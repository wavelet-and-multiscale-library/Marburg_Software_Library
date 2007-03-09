// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_A_EVALUATE_H
#define _WAVELETTL_A_EVALUATE_H

#include <utils/array1d.h>
#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <interval/a_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;
using MathTL::Array1D;

namespace WaveletTL
{
  /*!
    Evaluate a single primal generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].
  */
  template <unsigned int n>
  SampledMapping<1> evaluate(const ABasis<n>& basis,
                             const typename ABasis<n>::Index& lambda,
                             const int resolution);
  
  /*!
    Evaluate an arbitrary linear combination of primal
    wavelets on a dyadic subgrid of [0,1].
  */
  template <unsigned int n>
  SampledMapping<1> evaluate(const ABasis<n>& basis,
                             const InfiniteVector<double, typename ABasis<n>::Index>& coeffs,
                             const int resolution);

  /*!
    dummy routine
  */
  template <unsigned int n>
  SampledMapping<1> evaluate(const ABasis<n>& basis,
                             const InfiniteVector<double, typename ABasis<n>::Index>& coeffs,
                             const bool primal,
                             const int resolution)
  {
    return evaluate(basis, coeffs, resolution);
  }

  /*!
    point evaluation of (derivatives) of a single primal [A] generator
    or wavelet \psi_\lambda
  */
  template <unsigned int n>
  double evaluate(const ABasis<n>& basis, const unsigned int derivative,
                  const typename ABasis<n>::Index& lambda,
                  const double x);

  /*!
    point evaluation of (derivatives) of a single primal [A] generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <unsigned int n>
  void evaluate(const ABasis<n>& basis, const unsigned int derivative,
                const typename ABasis<n>::Index& lambda,
                const Array1D<double>& points, Array1D<double>& values);

  /*!
    point evaluation of 0-th and first derivative of a single primal [A] generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <unsigned int n>
  void evaluate(const ABasis<n>& basis,
                const typename ABasis<n>::Index& lambda,
                const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);

}

#include <interval/a_evaluate.cpp>

#endif // _WAVELETTL_A_EVALUATE_H

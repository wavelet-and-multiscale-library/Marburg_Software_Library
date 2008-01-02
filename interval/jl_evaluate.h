// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_JL_EVALUATE_H
#define _WAVELETTL_JL_EVALUATE_H

#include <utils/array1d.h>
#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <interval/jl_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;
using MathTL::Array1D;

namespace WaveletTL
{
  /*!
    Evaluate a single primal generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].
  */
  SampledMapping<1> evaluate(const JLBasis& basis,
			     const JLBasis::Index& lambda,
			     const int resolution);
  
  /*!
    Evaluate an arbitrary linear combination of primal
    wavelets on a dyadic subgrid of [0,1].
  */
  SampledMapping<1> evaluate(const JLBasis& basis,
 			     const InfiniteVector<double, JLBasis::Index>& coeffs,
 			     const int resolution);

  /*!
    dummy routine
  */
  SampledMapping<1> evaluate(const JLBasis& basis,
			     const InfiniteVector<double, JLBasis::Index>& coeffs,
			     const bool primal,
			     const int resolution)
  {
    return evaluate(basis, coeffs, resolution);
  }

  /*!
    point evaluation of (derivatives) of a single primal [JL] generator
    or wavelet \psi_\lambda
  */
  double evaluate(const JLBasis& basis, const unsigned int derivative,
 		  const JLBasis::Index& lambda,
 		  const double x);

  /*!
    point evaluation of (derivatives) of a single primal [JL] generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  void evaluate(const JLBasis& basis, const unsigned int derivative,
		const JLBasis::Index& lambda,
		const Array1D<double>& points, Array1D<double>& values);

  /*!
    point evaluation of 0-th and first derivative of a single primal [JL] generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  void evaluate(const JLBasis& basis,
		const JLBasis::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);
}

#include <interval/jl_evaluate.cpp>

#endif

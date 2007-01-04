// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DS_EVALUATE_H
#define _WAVELETTL_DS_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <interval/ds_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT, DSBiorthogonalizationMethod BIO> class DSBasis;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].

    We redirect this to a member function of basis.
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  SampledMapping<1> evaluate(const DSBasis<d,dT,BIO>& basis,
			     const typename DSBasis<d,dT,BIO>::Index& lambda,
			     const bool primal,
			     const int resolution)
  {
    return basis.evaluate(lambda, primal, resolution);
  }
  
  /*!
    Evaluate an arbitrary linear combination of primal or dual
    wavelets on a dyadic subgrid of [0,1].

    We redirect this to a member function of basis.
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  SampledMapping<1> evaluate(const DSBasis<d,dT,BIO>& basis,
			     const InfiniteVector<double, typename DSBasis<d,dT,BIO>::Index>& coeffs,
			     const bool primal,
			     const int resolution)
  {
    return basis.evaluate(coeffs, primal, resolution);
  }

  /*!
    point evaluation of (derivatives) of a single primal DKU generator
    or wavelet \psi_\lambda
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  double evaluate(const DSBasis<d,dT,BIO>& basis, const unsigned int derivative,
		  const typename DSBasis<d,dT,BIO>::Index& lambda,
		  const double x);

  /*!
    point evaluation of (derivatives) of a single primal DKU generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  void evaluate(const DSBasis<d,dT,BIO>& basis, const unsigned int derivative,
		const typename DSBasis<d,dT,BIO>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& values);

  /*!
    point evaluation of 0-th and first derivative of a single primal DKU generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  void evaluate(const DSBasis<d,dT,BIO>& basis,
		const typename DSBasis<d,dT,BIO>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);
}

#include <interval/ds_evaluate.cpp>

#endif

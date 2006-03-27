// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_P_EVALUATE_H
#define _WAVELETTL_P_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <interval/p_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT> class PBasis;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1].
  */
  template <int d, int dT>
  SampledMapping<1> evaluate(const PBasis<d,dT>& basis,
			     const typename PBasis<d,dT>::Index& lambda,
			     const bool primal,
			     const int resolution);
  
  /*!
    Evaluate an arbitrary linear combination of primal or dual
    wavelets on a dyadic subgrid of [0,1].
  */
  template <int d, int dT>
  SampledMapping<1> evaluate(const PBasis<d,dT>& basis,
			     const InfiniteVector<double, typename PBasis<d,dT>::Index>& coeffs,
			     const bool primal,
			     const int resolution);

  /*!
    point evaluation of (derivatives) of a single primal [P] generator
    or wavelet \psi_\lambda
  */
  template <int d, int dT>
  double evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
		  const typename PBasis<d,dT>::Index& lambda,
		  const double x);

  /*!
    point evaluation of (derivatives) of a single primal [P] generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <int d, int dT>
  void evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
		const typename PBasis<d,dT>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& values);

  /*!
    point evaluation of 0-th and first derivative of a single primal [P] generator
    or wavelet \psi_\lambda at several points simultaneously
  */
  template <int d, int dT>
  void evaluate(const PBasis<d,dT>& basis,
		const typename PBasis<d,dT>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);
}

#include <interval/p_evaluate.cpp>

#endif

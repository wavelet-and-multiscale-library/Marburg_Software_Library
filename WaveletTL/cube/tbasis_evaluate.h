// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TBASIS_EVALUATE_H
#define	_WAVELETTL_TBASIS_EVALUATE_H


#include <algebra/infinite_vector.h>
#include <cube/tbasis.h>
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
  template <class IBASIS, unsigned int DIM> class TensorBasis;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorBasis<IBASIS,DIM>& basis,
			       const typename TensorBasis<IBASIS,DIM>::Index& lambda,
			       const bool primal,
			       const int resolution);
  
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorBasis<IBASIS,DIM>& basis,
			       const int& lambdanum,
			       const bool primal,
			       const int resolution);

  /*!
    Evaluate an arbitrary linear combination of primal/dual wavelets
    on a dyadic subgrid of [0,1]^d.
  */
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorBasis<IBASIS,DIM>& basis,
			       const InfiniteVector<double, typename TensorBasis<IBASIS,DIM>::Index>& coeffs,
			       const bool primal,
			       const int resolution);
  
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const TensorBasis<IBASIS,DIM>& basis,
			       const InfiniteVector<double, int>& coeffs,
			       const bool primal,
			       const int resolution);
}

#include <cube/tbasis_evaluate.cpp>

#endif	/* _WAVELETTL_TBASIS_EVALUATE_H */


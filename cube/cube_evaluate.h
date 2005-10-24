// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_EVALUATE_H
#define _WAVELETTL_CUBE_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <cube/cube_basis.h>

using MathTL::SampledMapping;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM> class CubeBasis;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of [0,1]^2.
  */
  template <class IBASIS, unsigned int DIM>
  SampledMapping<DIM> evaluate(const CubeBasis<IBASIS,DIM>& basis,
			       const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			       const bool primal,
			       const int resolution);
}

#include <cube/cube_evaluate.cpp>

#endif

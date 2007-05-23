// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_EXPANSION_H
#define _WAVELETTL_CUBE_EXPANSION_H

#include <algebra/infinite_vector.h>
#include <cube/cube_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    For a given function, compute all integrals w.r.t. the primal
    or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
    - When integrating against the primal functions, the integrand has to be smooth
      to be accurately reproduced by the dual basis.
    - When integration against dual functions is specified,
      we integrate against the primal ones instead and multiply the resulting
      coefficients with the inverse of the primal gramian.

    Maybe a thresholding of the returned coefficients is helpful (e.g. for
    expansions of spline functions).
  */
  template <class IBASIS, unsigned int DIM>
  void
  expand(const Function<DIM>* f,
	 const CubeBasis<IBASIS,DIM>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index>& coeffs);
}

#include <cube/cube_expansion.cpp>

#endif

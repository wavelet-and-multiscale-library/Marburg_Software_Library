// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_P_EXPANSION_H
#define _WAVELETTL_P_EXPANSION_H

#include <algebra/infinite_vector.h>
#include <interval/p_basis.h>
#include <interval/spline_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT> class PBasis;

  /*!
    helper function, integrate a smooth function f against a
    primal [P] generator or wavelet
  */
  template <int d, int dT>
  double integrate(const Function<1>* f,
		   const PBasis<d,dT>& basis,
		   const typename PBasis<d,dT>::Index& lambda);

  /*!
    helper function, integrate two primal [P] generators or wavelets
    against each other (for the Gramian)
  */
  template <int d, int dT>
  double integrate(const PBasis<d,dT>& basis,
		   const typename PBasis<d,dT>::Index& lambda,
		   const typename PBasis<d,dT>::Index& mu);
  
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
  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const PBasis<d,dT>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename PBasis<d,dT>::Index>& coeffs);
}

#include <interval/p_expansion.cpp>

#endif

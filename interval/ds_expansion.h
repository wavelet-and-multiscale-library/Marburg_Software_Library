// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DS_EXPANSION_H
#define _WAVELETTL_DS_EXPANSION_H

#include <algebra/infinite_vector.h>
#include <interval/ds_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT, DSBiorthogonalizationMethod BIO> class DSBasis;

  /*!
    helper function, integrate a smooth function f against a
    primal DS generator or wavelet
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  double integrate(const Function<1>* f,
		   const DSBasis<d,dT,BIO>& basis,
		   const typename DSBasis<d,dT,BIO>::Index& lambda);

  /*!
    helper function, integrate two primal DS generators or wavelets
    against each other (for the Gramian)
  */
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  double integrate(const DSBasis<d,dT,BIO>& basis,
		   const typename DSBasis<d,dT,BIO>::Index& lambda,
		   const typename DSBasis<d,dT,BIO>::Index& mu);
  
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
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  void
  expand(const Function<1>* f,
	 const DSBasis<d,dT,BIO>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename DSBasis<d,dT,BIO>::Index>& coeffs);
}

#include <interval/ds_expansion.cpp>

#endif

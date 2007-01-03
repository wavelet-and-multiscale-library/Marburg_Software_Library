// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_EXPANSION_H
#define _WAVELETTL_SPLINE_EXPANSION_H

#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor> class SplineBasis;
  template <int d, int dT> class SplineBasis<d,dT,DS_construction>;

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
  template <int d, int dT, SplineBasisFlavor flavor>
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,flavor>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs);

  /*!
    analogous routine for Vector<double> output
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  void
  expand(const Function<1>* f,
	 const SplineBasis<d,dT,flavor>& basis,
	 const bool primal,
	 const int jmax,
	 Vector<double>& coeffs);
}

#include <interval/spline_expansion.cpp>

#endif

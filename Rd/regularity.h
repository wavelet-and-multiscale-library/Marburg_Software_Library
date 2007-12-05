// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_VILLEMOES_H
#define _WAVELETTL_VILLEMOES_H

#include <iostream>
#include <algebra/multi_laurent_polynomial.h>
#include <utils/multiindex.h>

using MathTL::MultivariateLaurentPolynomial;
using MathTL::MultiIndex;

namespace WaveletTL
{
  /*!
    For a given univariate mask a, compute the autocorrelation mask

      b(k) = \sum_m a(k+m)*a(m)/2

    The autocorrelation mask is needed for the calculation of the Sobolev smoothness
    of a given refinable function.
   */
  template <class MASK>
  class AutocorrelationMask
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    AutocorrelationMask();
  };

  /*!
    Compute the critical Sobolev regularity of a given univariate
    refinable function, see [V] or [JZ] for details.
    The template parameter class MASK should provide a method Strang_Fix_order().

    UNDER CONSTRUCTION!!!

    References:

    [JZ] Jia/Zhang: Spectral properties of the transition operator associated to a
         multivariate refinement equation,
	 Linear Algebra and its applications 292 (1999), 155-178
    [V]  Villemoes: Wavelet analysis of two-scale refinement equations,
         SIAM J. Math. Anal. 25 (1994), 1433-1460
  */
  template <class MASK>
  double
  Sobolev_regularity();
}

#include <Rd/regularity.cpp>

#endif

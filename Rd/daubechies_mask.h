// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DAUBECHIES_MASK_H
#define _WAVELETTL_DAUBECHIES_MASK_H

#include <iostream>
#include <algebra/multi_laurent_polynomial.h>
#include <utils/multiindex.h>

using MathTL::MultivariateLaurentPolynomial;
using MathTL::MultiIndex;

namespace WaveletTL
{
  /*!
    The masks of the orthogonal refinable functions of order N constructed
    by I. Daubechies in [D], N being the polynomial exactness of the generators
    as well as the number of vanishing moments of the wavelets.

    [D]  Daubechies: Orthonormal bases of compactly supported wavelets,
         Commun. Pure Appl. Math. 41(1988), 909--996
    [SY] Shann/Yen: Exact solutions for Daubechies orthonormal scaling coefficients,
         Technical Report, 1997
    [SY] Shann/Yen: On the exact values of orthonormal scaling coefficients of lengths 8 and 10,
         Appl. Comput. Harmon. Anal. 6(1999), 109--112
  */
  template <int N>
  class DaubechiesMask
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    DaubechiesMask();
  };
}

#include <Rd/daubechies_mask.cpp>

#endif

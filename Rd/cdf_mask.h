// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_MASK_H
#define _WAVELETTL_CDF_MASK_H

#include <iostream>
#include <algebra/multi_laurent_polynomial.h>

using MathTL::MultivariateLaurentPolynomial;

namespace WaveletTL
{
  /*!
    primal mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d>
  class CDFMask_primal
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    CDFMask_primal();
  };

  /*!
    dual mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d, int dt>
  class CDFMask_dual
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    CDFMask_dual();
  };
}

// include implementation
#include <Rd/cdf_mask.cpp>

#endif

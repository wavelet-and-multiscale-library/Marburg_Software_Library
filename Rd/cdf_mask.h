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
#include <algebra/laurent_polynomial.h>

using MathTL::LaurentPolynomial;

namespace WaveletTL
{
  /*!
    primal mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <unsigned int d>
  class CDFMask_primal
    : public virtual LaurentPolynomial<double>
  {
  public:
    CDFMask_primal();
  };

  /*!
    dual mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <unsigned int d, unsigned int dt>
  class CDFMask_dual
    : public virtual LaurentPolynomial<double>
  {
  public:
    CDFMask_dual();
  };
}

// include implementation
#include <Rd/cdf_mask.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TRIVIAL_MASK_H
#define _WAVELETTL_TRIVIAL_MASK_H

#include <iostream>
#include <algebra/laurent_polynomial.h>

using MathTL::LaurentPolynomial;

namespace WaveletTL
{
  class TrivialMask
    : public virtual LaurentPolynomial<double>
  {
  public:
    TrivialMask()
      : LaurentPolynomial<double>(1)
    {
    }
  };
}

#endif

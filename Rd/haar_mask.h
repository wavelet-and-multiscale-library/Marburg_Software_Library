// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_HAAR_MASK_H
#define _WAVELETTL_HAAR_MASK_H

#include <iostream>
#include <algebra/laurent_polynomial.h>

using MathTL::LaurentPolynomial;

namespace WaveletTL
{
  class HaarMask
    : public virtual LaurentPolynomial<double>
    {
    public:
      HaarMask()
	{
	  set_coefficient(0, 1);
	  set_coefficient(1, 1);
	}
    };
}

#endif

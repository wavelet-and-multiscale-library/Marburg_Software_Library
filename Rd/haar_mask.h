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
#include <algebra/multi_laurent_polynomial.h>
#include <utils/multiindex.h>

using MathTL::MultivariateLaurentPolynomial;
using MathTL::MultiIndex;

namespace WaveletTL
{
  class HaarMask
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    HaarMask()
    {
      set_coefficient(MultiIndex<int, 1>(0), 1);
      set_coefficient(MultiIndex<int, 1>(1), 1);
    }

    /*!
      Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    */
    static unsigned int Strang_Fix_order() { return 1; }

    /*!
      critical Sobolev regularity
    */
    static double regularity() { return 0.5; }
  };
}

#endif

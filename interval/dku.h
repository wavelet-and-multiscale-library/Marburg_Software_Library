// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DKU_H
#define _WAVELETTL_DKU_H

#include <iostream>

namespace WaveletTL
{
  /*!
    Template class for the wavelet bases on the interval introduced in [DKU].

    References:
    [DKU]
  */
  template <unsigned int d, unsigned int dt>
  class DKUBasis
  {
  public:
  };
}

#include <interval/dku.cpp>

#endif

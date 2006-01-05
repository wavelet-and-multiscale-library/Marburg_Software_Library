// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_EVALUATE_H
#define _WAVELETTL_CDF_EVALUATE_H

namespace WaveletTL
{
  /*!
    point evaluation of primal CDF generators/wavelets
  */
  template <int d, int dT>
  double evaluate(const unsigned int derivative,
		  const RIndex& lambda,
		  const double x);
}

#include <Rd/cdf_evaluate.cpp>

#endif

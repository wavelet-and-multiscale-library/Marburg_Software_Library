// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DKU_EVALUATE_H
#define _WAVELETTL_DKU_EVALUATE_H

#include <interval/dku_basis.h>

namespace WaveletTL
{
  /*!
    point evaluation of (derivatives) of a single primal DKU generator
    or wavelet \psi_\lambda
   */
  template <int d, int dT>
  double evaluate(const DKUBasis<d, dT>& basis, const unsigned int derivative,
		  const typename DKUBasis<d, dT>::Index& lambda,
		  const double x);
}

#include <interval/dku_evaluate.cpp>

#endif

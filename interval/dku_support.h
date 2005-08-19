// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DKU_SUPPORT_H
#define _WAVELETTL_DKU_SUPPORT_H

namespace WaveletTL
{
  /*!
    compute an interval 2^{-j}[k1,k2] which contains the support of a
    single primal DKU generator or wavelet \psi_\lambda
    (where j=lambda.j()+lambda.e())
   */
  template <int d, int dT>
  void support(const DKUBasis<d, dT>& basis,
	       const typename DKUBasis<d, dT>::Index& lambda,
	       int& k1, int& k2);
}

#include <interval/dku_support.cpp>

#endif

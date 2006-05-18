// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_SUPPORT_H
#define _WAVELETTL_LDOMAIN_SUPPORT_H

#include <list>
#include <set>

namespace WaveletTL
{
  template <class IBASIS> class LDomainBasis;

  /*!
    Compute a set which contains the support of a single primal generator
    or wavelet psi_lambda.
  */
  template <class IBASIS>
  void support(const LDomainBasis<IBASIS>& basis,
	       const typename LDomainBasis<IBASIS>::Index& lambda,
	       typename LDomainBasis<IBASIS>::Support& supp);
}

#include <Ldomain/ldomain_support.cpp>

#endif

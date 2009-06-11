// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_SUPPORT_H
#define _WAVELETTL_LDOMAIN_SUPPORT_H

#include <list>
#include <set>
#include <algebra/infinite_vector.h>

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

  /*!
    Compute the support intersection of a wavelet psi_lambda with the
    support of another wavelet psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IBASIS>
  bool intersect_supports(const LDomainBasis<IBASIS>& basis,
			  const typename LDomainBasis<IBASIS>::Index& lambda,
			  const typename LDomainBasis<IBASIS>::Support& supp_mu,
			  typename LDomainBasis<IBASIS>::Support& supp);

  /*!
    Compute the support intersection of two wavelets psi_lambda, psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IBASIS>
  bool intersect_supports(const LDomainBasis<IBASIS>& basis,
			  const typename LDomainBasis<IBASIS>::Index& lambda,
			  const typename LDomainBasis<IBASIS>::Index& mu,
			  typename LDomainBasis<IBASIS>::Support& supp);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <class IBASIS>
  void intersecting_wavelets(const LDomainBasis<IBASIS>& basis,
			     const typename LDomainBasis<IBASIS>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename LDomainBasis<IBASIS>::Index>& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
  */
  template <class IBASIS>
  bool intersect_singular_support(const LDomainBasis<IBASIS>& basis,
				  const typename LDomainBasis<IBASIS>::Index& lambda,
				  const typename LDomainBasis<IBASIS>::Index& nu);

}

#include <Ldomain/ldomain_support.cpp>

#endif

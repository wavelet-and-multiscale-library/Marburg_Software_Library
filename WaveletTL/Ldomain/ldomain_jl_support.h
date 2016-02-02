// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_JL_SUPPORT_H
#define _WAVELETTL_LDOMAIN_JL_SUPPORT_H

#include <list>
#include <set>
#include <algebra/infinite_vector.h>
#include <Ldomain/ldomain_jl_basis.h>

namespace WaveletTL
{
  typedef LDomainJLBasis::Index Index;
  typedef LDomainJLBasis::Support Support;
  
  /*!
    Compute a set which contains the support of a single primal generator
    or wavelet psi_lambda.
    To be fast, we always compute a square that is not (!) necessarily a subset
    of Omega. For the interior generators/wavelets, the support is correct.
    Only for the boundary generators/wavelets, we have to be careful in
    point evaluations.
  */
  void support(const LDomainJLBasis& basis,
 	       const Index& lambda,
 	       Support& supp);
  
  /*!
    Compute the support intersection of a wavelet psi_lambda with the
    support of another wavelet psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  bool intersect_supports(const LDomainJLBasis& basis,
			  const Index& lambda,
			  const Support& supp_mu,
			  Support& supp);

  /*!
    Compute the support intersection of two wavelets psi_lambda, psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  bool intersect_supports(const LDomainJLBasis& basis,
			  const Index& lambda,
			  const Index& mu,
			  Support& supp);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  void intersecting_wavelets(const LDomainJLBasis& basis,
			     const Index& lambda,
			     const int j, const bool generators,
			     std::list<Index>& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
  */
  bool intersect_singular_support(const LDomainJLBasis& basis,
				  const Index& lambda,
				  const Index& nu);
  
}

#include <Ldomain/ldomain_jl_support.cpp>

#endif

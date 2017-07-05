// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_FRAME_SUPPORT_H
#define _WAVELETTL_LDOMAIN_FRAME_SUPPORT_H

#include <list>
#include <set>
#include <algebra/infinite_vector.h>

namespace WaveletTL
{
  template <class IFRAME> class LDomainFrame;

  /*!
    Compute a set which contains the support of a single primal generator
    or quarklet psi_lambda.
  */
  template <class IFRAME>
  void support(const LDomainFrame<IFRAME>& frame,
	       const typename LDomainFrame<IFRAME>::Index& lambda,
	       typename LDomainFrame<IFRAME>::Support& supp);

  /*!
    Compute the support intersection of a quarklet psi_lambda with the
    support of another quarklet psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IFRAME>
  bool intersect_supports(const LDomainFrame<IFRAME>& frame,
			  const typename LDomainFrame<IFRAME>::Index& lambda,
			  const typename LDomainFrame<IFRAME>::Support& supp_mu,
			  typename LDomainFrame<IFRAME>::Support& supp);

  /*!
    Compute the support intersection of two quarklets psi_lambda, psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IFRAME>
  bool intersect_supports(const LDomainFrame<IFRAME>& frame,
			  const typename LDomainFrame<IFRAME>::Index& lambda,
			  const typename LDomainFrame<IFRAME>::Index& mu,
			  typename LDomainFrame<IFRAME>::Support& supp);

  /*!
    For a given quarklet \psi_\lambda, compute all generators/quarklets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <class IFRAME>
  void intersecting_quarklets(const LDomainFrame<IFRAME>& frame,
			     const typename LDomainFrame<IFRAME>::Index& lambda,
			     const typename LDomainFrameIndex<IFRAME>::level_type& j,
			     std::list<typename LDomainFrame<IFRAME>::Index>& intersecting,
                             const typename LDomainFrameIndex<IFRAME>::polynomial_type& p);

  /*!
    Decide whether the support of a given (primal) generator/quarklet \psi_\lambda
    intersects the singular support of another (primal) generator/quarklet \psi_\nu.
  */
  template <class IFRAME>
  bool intersect_singular_support(const LDomainFrame<IFRAME>& frame,
				  const typename LDomainFrame<IFRAME>::Index& lambda,
				  const typename LDomainFrame<IFRAME>::Index& nu);

}

#include <Ldomain/ldomain_frame_support.cpp>

#endif

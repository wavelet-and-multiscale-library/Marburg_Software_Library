// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DOMAIN_FRAME_SUPPORT_H
#define _WAVELETTL_DOMAIN_FRAME_SUPPORT_H

#include <list>
#include <set>
#include <algebra/infinite_vector.h>
#include <general_domain/domain_frame_index.h>

namespace WaveletTL
{
  template <class IFRAME, int NPATCHES> class DomainFrame;

  /*!
    Compute a set which contains the support of a single primal generator
    or quarklet psi_lambda.
  */
  template <class IFRAME, int NPATCHES>
  void support(const DomainFrame<IFRAME, NPATCHES>& frame,
	       const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
	       typename DomainFrame<IFRAME, NPATCHES>::Support& supp);
  
  template <class IFRAME, int NPATCHES>
  void support(const DomainFrame<IFRAME, NPATCHES>& frame,
	       const int& lambda_num,
	       typename DomainFrame<IFRAME, NPATCHES>::Support& supp);

  template <class IFRAME, int NPATCHES>
  bool intersect_supports(const DomainFrame<IFRAME, NPATCHES>& frame,
			  const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
			  const typename DomainFrame<IFRAME, NPATCHES>::Index& mu);
  
  
  
  /*!
    Compute the support intersection of a quarklet psi_lambda with the
    support of another quarklet psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IFRAME, int NPATCHES>
  bool intersect_supports(const DomainFrame<IFRAME, NPATCHES>& frame,
			  const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
			  const typename DomainFrame<IFRAME, NPATCHES>::Support& supp_mu,
			  typename DomainFrame<IFRAME, NPATCHES>::Support& supp);
  
  template <class IFRAME, int NPATCHES>
  bool intersect_supports(const DomainFrame<IFRAME, NPATCHES>& frame,
			  const int& lambda_num,
			  const typename DomainFrame<IFRAME, NPATCHES>::Support& supp_mu,
			  typename DomainFrame<IFRAME, NPATCHES>::Support& supp);

  /*!
    Compute the support intersection of two quarklets psi_lambda, psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IFRAME, int NPATCHES>
  bool intersect_supports(const DomainFrame<IFRAME, NPATCHES>& frame,
			  const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
			  const typename DomainFrame<IFRAME, NPATCHES>::Index& mu,
			  typename DomainFrame<IFRAME, NPATCHES>::Support& supp);
  
  template <class IFRAME, int NPATCHES>
  bool intersect_supports(const DomainFrame<IFRAME, NPATCHES>& frame,
			  const int& lambda_num,
			  const int& mu_num,
			  typename DomainFrame<IFRAME, NPATCHES>::Support& supp);

  /*!
    For a given quarklet \psi_\lambda, compute all generators/quarklets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <class IFRAME, int NPATCHES>
  void intersecting_quarklets(const DomainFrame<IFRAME, NPATCHES>& frame,
			     const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
			     const typename DomainFrameIndex<IFRAME, NPATCHES>::level_type& j,
                             std::list<int>& intersecting,
                             const typename DomainFrameIndex<IFRAME, NPATCHES>::polynomial_type& p);
  
  template <class IFRAME, int NPATCHES>
  void intersecting_quarklets(const DomainFrame<IFRAME, NPATCHES>& frame,
			     const int& lambda_num,
			     const typename DomainFrameIndex<IFRAME, NPATCHES>::level_type& j,
			     std::list<int>& intersecting,
                             const typename DomainFrameIndex<IFRAME, NPATCHES>::polynomial_type& p);
  
  

  /*!
    Decide whether the support of a given (primal) generator/quarklet \psi_\lambda
    intersects the singular support of another (primal) generator/quarklet \psi_\nu.
  */
  template <class IFRAME, int NPATCHES>
  bool intersect_singular_support(const DomainFrame<IFRAME, NPATCHES>& frame,
				  const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
				  const typename DomainFrame<IFRAME, NPATCHES>::Index& nu);

}

#include <general_domain/domain_frame_support.cpp>

#endif

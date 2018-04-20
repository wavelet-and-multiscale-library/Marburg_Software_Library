// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2018                                            |
// | Thorsten Philipp Keding, Alexander Sieber                          |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SLITDOMAIN_FRAME_SUPPORT_H
#define _WAVELETTL_SLITDOMAIN_FRAME_SUPPORT_H

#include <list>
#include <set>
#include <algebra/infinite_vector.h>

namespace WaveletTL
{
  template <class IFRAME> class SlitDomainFrame;

  /*!
    Compute a set which contains the support of a single primal generator
    or quarklet psi_lambda.
  */
  template <class IFRAME>
  void support(const SlitDomainFrame<IFRAME>& frame,
	       const typename SlitDomainFrame<IFRAME>::Index& lambda,
	       typename SlitDomainFrame<IFRAME>::Support& supp);
  
  template <class IFRAME>
  void support(const SlitDomainFrame<IFRAME>& frame,
	       const int& lambda_num,
	       typename SlitDomainFrame<IFRAME>::Support& supp);

  template <class IFRAME>
  bool intersect_supports(const SlitDomainFrame<IFRAME>& frame,
			  const typename SlitDomainFrame<IFRAME>::Index& lambda,
			  const typename SlitDomainFrame<IFRAME>::Index& mu);
  
  
  
  /*!
    Compute the support intersection of a quarklet psi_lambda with the
    support of another quarklet psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IFRAME>
  bool intersect_supports(const SlitDomainFrame<IFRAME>& frame,
			  const typename SlitDomainFrame<IFRAME>::Index& lambda,
			  const typename SlitDomainFrame<IFRAME>::Support& supp_mu,
			  typename SlitDomainFrame<IFRAME>::Support& supp);
  
  template <class IFRAME>
  bool intersect_supports(const SlitDomainFrame<IFRAME>& frame,
			  const int& lambda_num,
			  const typename SlitDomainFrame<IFRAME>::Support& supp_mu,
			  typename SlitDomainFrame<IFRAME>::Support& supp);

  /*!
    Compute the support intersection of two quarklets psi_lambda, psi_mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IFRAME>
  bool intersect_supports(const SlitDomainFrame<IFRAME>& frame,
			  const typename SlitDomainFrame<IFRAME>::Index& lambda,
			  const typename SlitDomainFrame<IFRAME>::Index& mu,
			  typename SlitDomainFrame<IFRAME>::Support& supp);
  
  template <class IFRAME>
  bool intersect_supports(const SlitDomainFrame<IFRAME>& frame,
			  const int& lambda_num,
			  const int& mu_num,
			  typename SlitDomainFrame<IFRAME>::Support& supp);

  /*!
    For a given quarklet \psi_\lambda, compute all generators/quarklets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <class IFRAME>
  void intersecting_quarklets(const SlitDomainFrame<IFRAME>& frame,
			     const typename SlitDomainFrame<IFRAME>::Index& lambda,
			     const typename SlitDomainFrameIndex<IFRAME>::level_type& j,
                             std::list<int>& intersecting,
                             const typename SlitDomainFrameIndex<IFRAME>::polynomial_type& p);
  
  template <class IFRAME>
  void intersecting_quarklets(const SlitDomainFrame<IFRAME>& frame,
			     const int& lambda_num,
			     const typename SlitDomainFrameIndex<IFRAME>::level_type& j,
			     std::list<int>& intersecting,
                             const typename SlitDomainFrameIndex<IFRAME>::polynomial_type& p);
  
  

  /*!
    Decide whether the support of a given (primal) generator/quarklet \psi_\lambda
    intersects the singular support of another (primal) generator/quarklet \psi_\nu.
  */
  template <class IFRAME>
  bool intersect_singular_support(const SlitDomainFrame<IFRAME>& frame,
				  const typename SlitDomainFrame<IFRAME>::Index& lambda,
				  const typename SlitDomainFrame<IFRAME>::Index& nu);

}

#include <slitdomain/slitdomain_frame_support.cpp>

#endif


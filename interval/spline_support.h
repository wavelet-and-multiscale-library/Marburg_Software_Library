// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SPLINE_SUPPORT_H
#define _WAVELETTL_SPLINE_SUPPORT_H

#include <list>
#include <set>
#include <interval/spline_basis_flavor.h>

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor> class SplineBasis;
  template <int d, int dT> class SplineBasis<d,dT,DS_construction>;

  /*!
    Compute an interval 2^{-j}[k1,k2] which contains the support of a
    single primal [P] generator or wavelet \psi_\lambda.
    (j == lambda.j()+lambda.e() is neglected for performance reasons)
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  void support(const SplineBasis<d,dT,flavor>& basis,
	       const typename SplineBasis<d,dT,flavor>::Index& lambda,
	       int& k1, int& k2);
  
  /*!
    Decide whether the support of a generator/wavelet \psi_\lambda
    has a nontrivial intersection with a given dyadic interval 2^{-m}[a,b]
    and compute it in the form 2^{-j}[k1,k2]. If the return value is false,
    k1 and k2 will have no meaningful value for performance reasons.
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  bool intersect_supports(const SplineBasis<d,dT,flavor>& basis,
			  const typename SplineBasis<d,dT,flavor>::Index& lambda,
			  const int m, const int a, const int b,
			  int& j, int& k1, int& k2);
  
  /*!
    Decide whether the supports of two generators/wavelets \psi_\lambda and
    \psi_\nu have an intersection of positive measure and compute it
    in the form 2^{-j}[k1,k2]. If the return value is false, the computed
    support intersection will have no meaningful values, for performance reasons.
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  bool intersect_supports(const SplineBasis<d,dT,flavor>& basis,
			  const typename SplineBasis<d,dT,flavor>::Index& lambda,
			  const typename SplineBasis<d,dT,flavor>::Index& nu,
			  typename SplineBasis<d,dT,flavor>::Support& supp);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  void intersecting_wavelets(const SplineBasis<d,dT,flavor>& basis,
			     const typename SplineBasis<d,dT,flavor>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename SplineBasis<d,dT,flavor>::Index>& intersecting);
  
  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection;
    return those intersections
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  void intersecting_wavelets
  (const SplineBasis<d,dT,flavor>& basis,
   const typename SplineBasis<d,dT,flavor>::Index& lambda,
   const int j, const bool generators,
   std::list<std::pair<typename SplineBasis<d,dT,flavor>::Index, typename SplineBasis<d,dT,flavor>::Support> >& intersecting);
  
  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of a union 2^{-m}[a,b] of dyadic unit intervals.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
    The routine should only be called in the case lambda.j()+lambda.e() >= m.
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  bool intersect_singular_support(const SplineBasis<d,dT,flavor>& basis,
				  const typename SplineBasis<d,dT,flavor>::Index& lambda,
				  const int m, const int a, const int b,
				  int& j, int& k1, int& k2);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  bool intersect_singular_support(const SplineBasis<d,dT,flavor>& basis,
				  const typename SplineBasis<d,dT,flavor>::Index& lambda,
				  const typename SplineBasis<d,dT,flavor>::Index& nu);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
  */
  template <int d, int dT, SplineBasisFlavor flavor>
  bool intersect_singular_support(const SplineBasis<d,dT,flavor>& basis,
				  const typename SplineBasis<d,dT,flavor>::Index& lambda,
				  const typename SplineBasis<d,dT,flavor>::Index& nu,
				  int& j, int& k1, int& k2);

}

#include <interval/spline_support.cpp>

#endif

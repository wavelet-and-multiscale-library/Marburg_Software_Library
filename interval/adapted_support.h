// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_ADAPTED_SUPPORT_H
#define _WAVELETTL_ADAPTED_SUPPORT_H

#include <interval/adapted_basis.h>

#include <list>
#include <set>

namespace WaveletTL
{
  template <class IBASIS> class AdaptedBasis;

  /*!
    Compute an interval 2^{-j}[k1,k2] which contains the support of a
    single primal [P] generator or wavelet \psi_\lambda.
    (j == lambda.j()+lambda.e() is neglected for performance reasons)
  */
  template <class IBASIS>
  void support(const AdaptedBasis<IBASIS>& basis,
               const typename AdaptedBasis<IBASIS>::Index& lambda,
               int& k1, int& k2);
  
  /*!
    Decide whether the support of a generator/wavelet \psi_\lambda
    has a nontrivial intersection with a given dyadic interval 2^{-m}[a,b]
    and compute it in the form 2^{-j}[k1,k2]. If the return value is false,
    k1 and k2 will have no meaningful value for performance reasons.
  */
  template <class IBASIS>
  bool intersect_supports(const AdaptedBasis<IBASIS>& basis,
                          const typename AdaptedBasis<IBASIS>::Index& lambda,
                          const int m, const int a, const int b,
                          int& j, int& k1, int& k2);
  
  /*!
    Decide whether the supports of two generators/wavelets \psi_\lambda and
    \psi_\nu have an intersection of positive measure and compute it
    in the form 2^{-j}[k1,k2]. If the return value is false, the computed
    support intersection will have no meaningful values, for performance reasons.
  */
  template <class IBASIS>
  bool intersect_supports(const AdaptedBasis<IBASIS>& basis,
                          const typename AdaptedBasis<IBASIS>::Index& lambda,
                          const typename AdaptedBasis<IBASIS>::Index& nu,
                          typename AdaptedBasis<IBASIS>::Support& supp);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <class IBASIS>
  void intersecting_wavelets(const AdaptedBasis<IBASIS>& basis,
                             const typename AdaptedBasis<IBASIS>::Index& lambda,
                             const int j, const bool generators,
                             std::list<typename AdaptedBasis<IBASIS>::Index>& intersecting);
  
  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection;
    return those intersections
  */
  template <class IBASIS>
  void intersecting_wavelets(const AdaptedBasis<IBASIS>& basis,
                             const typename AdaptedBasis<IBASIS>::Index& lambda,
                             const int j, const bool generators,
                             std::list<std::pair<typename AdaptedBasis<IBASIS>::Index, typename AdaptedBasis<IBASIS>::Support> >& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of a union 2^{-m}[a,b] of dyadic unit intervals.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
    The routine should only be called in the case lambda.j()+lambda.e() >= m.
  */
  template <class IBASIS>
  bool intersect_singular_support(const AdaptedBasis<IBASIS>& basis,
                                  const typename AdaptedBasis<IBASIS>::Index& lambda,
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
  template <class IBASIS>
  bool intersect_singular_support(const AdaptedBasis<IBASIS>& basis,
                                  const typename AdaptedBasis<IBASIS>::Index& lambda,
                                  const typename AdaptedBasis<IBASIS>::Index& nu);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
  */
  template <class IBASIS>
  bool intersect_singular_support(const AdaptedBasis<IBASIS>& basis,
                                  const typename AdaptedBasis<IBASIS>::Index& lambda,
                                  const typename AdaptedBasis<IBASIS>::Index& nu,
                                  int& j, int& k1, int& k2);

}

#include <interval/adapted_support.cpp>

#endif // _WAVELETTL_ADAPTED_SUPPORT_H

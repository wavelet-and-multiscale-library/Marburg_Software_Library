// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_S_SUPPORT_H
#define _WAVELETTL_S_SUPPORT_H

#include <list>
#include <set>

namespace WaveletTL
{
  /*!
    Compute an interval 2^{-j}[k1,k2] which contains the support of a
    single primal [S] generator or wavelet \psi_\lambda.
    Note that wavelets have granularity lambda.j()+1.
    \param basis underlying basis
    \param lambda index of the generator/wavelet to calculate the support of
    \param k1 left support interval boundary in terms of 2^{-j}
    \param k2 right support interval boundary in terms of 2^{-j}
  */
  inline
  void support(const SBasis& basis,
               const SBasis::Index& lambda,
               int& k1, int& k2);

  /*!
    Decide whether the support of a generator/wavelet \psi_\lambda
    has a nontrivial intersection with a given dyadic interval 2^{-m}[a,b]
    and compute it in the form 2^{-j}[k1,k2]. If the return value is false,
    k1 and k2 will have no meaningful value for performance reasons.
  */
  bool intersect_supports(const SBasis& basis,
                          const SBasis::Index& lambda,
                          const int m, const int a, const int b,
                          int& j, int& k1, int& k2);

  /*!
    Decide whether the supports of two generators/wavelets \psi_\lambda and
    \psi_\nu have an intersection of positive measure and compute it
    in the form 2^{-j}[k1,k2]. If the return value is false, the computed
    support intersection will have no meaningful values, for performance reasons.
  */
  inline
  bool intersect_supports(const SBasis& basis,
                          const SBasis::Index& lambda,
                          const SBasis::Index& nu,
                          SBasis::Support& supp);

  /*!
    For a given wavelet \psi_\lambda, compute translation indices k1, k2, such that
    the respective supports all generators/wavelets \psi_\nu with level |\nu|=j
    and translation index k1 <= k <= k2 have a nontrivial intersection    
    \param basis underlying basis
    \param lambda index of the generator/wavelet to calculate intersecting functions of
    \param j level on which intersecting functions are searched for
    \param generators decide whether to search for intersecting generators or wavelets
    \param k1 starting index (on the left side)
    \param k2 ending index (on the right side)
  */
  void intersecting_wavelets(const SBasis& basis,
                             const SBasis::Index& lambda,
                             const int j, const bool generators,
                             int& k1, int& k2);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
    \param basis underlying basis
    \param lambda index of the generator/wavelet to calculate intersecting functions of
    \param j level on which intersecting functions are searched for
    \param generators decide whether to search for intersecting generators or wavelets
    \param intersecting list to be filled with the indices of the intersecting functions
  */
  void intersecting_wavelets(const SBasis& basis,
                             const SBasis::Index& lambda,
                             const int j, const bool generators,
                             std::list<SBasis::Index>& intersecting);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection;
    return those intersections
  */
  void intersecting_wavelets
    (const SBasis& basis,
     const SBasis::Index& lambda,
     const int j, const bool generators,
     std::list<std::pair<SBasis::Index, SBasis::Support> >& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
  */
  bool intersect_singular_support(const SBasis& basis,
                                  const SBasis::Index& lambda,
                                  const SBasis::Index& nu);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
  */
  bool intersect_singular_support(const SBasis& basis,
                                  const SBasis::Index& lambda,
                                  const SBasis::Index& nu,
                                  int& j, int& k1, int& k2);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of a union 2^{-m}[a,b] of dyadic unit intervals.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
    The routine should only be called in the case lambda.j()+lambda.e() >= m.
  */
  bool intersect_singular_support(const SBasis& basis,
                                  const SBasis::Index& lambda,
                                  const int m, const int a, const int b,
                                  int& j, int& k1, int& k2);
}

#include <interval/s_support.cpp>

#endif // _WAVELETTL_S_SUPPORT_H

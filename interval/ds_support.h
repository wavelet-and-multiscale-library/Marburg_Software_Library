// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DS_SUPPORT_H
#define _WAVELETTL_DS_SUPPORT_H

#include <list>
#include <set>

namespace WaveletTL
{
  template <int d, int dT> class DSBasis;

  /*!
    Compute an interval 2^{-j}[k1,k2] which contains the support of a
    single primal DS generator or wavelet \psi_\lambda.
    (j == lambda.j()+lambda.e() is neglected for performance reasons)
   */
  template <int d, int dT>
  void support(const DSBasis<d,dT>& basis,
	       const typename DSBasis<d,dT>::Index& lambda,
	       int& k1, int& k2);

  /*!
    Decide whether the support of a generator/wavelet \psi_\lambda
    has a nontrivial intersection with a given dyadic interval 2^{-m}[a,b]
    and compute it in the form 2^{-j}[k1,k2]. If the return value is false,
    k1 and k2 will have no meaningful value for performance reasons.
  */
  template <int d, int dT>
  bool intersect_supports(const DSBasis<d,dT>& basis,
			  const typename DSBasis<d,dT>::Index& lambda,
			  const int m, const int a, const int b,
			  int& j, int& k1, int& k2);
  
  /*!
    Decide whether the supports of two generators/wavelets \psi_\lambda and
    \psi_\nu have an intersection of positive measure and compute it
    in the form 2^{-j}[k1,k2]. If the return value is false, k1 and k2 will
    have no meaningful value for performance reasons.
  */
  template <int d, int dT>
  bool intersect_supports(const DSBasis<d,dT>& basis,
			  const typename DSBasis<d,dT>::Index& lambda,
			  const typename DSBasis<d,dT>::Index& nu,
			  int& j, int& k1, int& k2);

  /*!
    a type for the support intervals 2^{-j}[k1,k2] we need here
  */
  typedef struct {
    int j;
    int k1;
    int k2;
  } Support1D;

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection;
    return those intersections
    (on the coarsest level j==j0, both generators and wavelets are used,
    whereas for j>j0, there are only wavelets)
  */
  template <int d, int dT>
  void intersecting_wavelets(const DSBasis<d,dT>& basis,
			     const typename DSBasis<d,dT>::Index& lambda,
			     const int j,
			     std::list<std::pair<typename DSBasis<d,dT>::Index, Support1D> >& intersecting);
}

#include <interval/ds_support.cpp>

#endif

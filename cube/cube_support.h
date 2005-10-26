// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_SUPPORT_H
#define _WAVELETTL_CUBE_SUPPORT_H

#include <list>
#include <set>

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM> class CubeBasis;

  /*!
    For a given interval basis IBASIS, compute a cube
    2^{-j}<a,b> = 2^{-j}[a_1,b_1]x...x[a_n,b_n]
    which contains the support of a single primal cube generator
    or wavelet psi_lambda.
  */
  template <class IBASIS, unsigned int DIM>
  void support(const CubeBasis<IBASIS,DIM>& basis,
	       const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	       typename CubeBasis<IBASIS,DIM>::Support& supp);
    
  /*!
    For a given interval basis, compute a cube
    2^{-j}<a,b> = 2^{-j}[a_1,b_1]x...x[a_n,b_n]
    representing the intersection of the support cubes
    corresponding to the indices lambda and mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IBASIS, unsigned int DIM>
  bool intersect_supports(const CubeBasis<IBASIS,DIM>& basis,
			  const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			  const typename CubeBasis<IBASIS,DIM>::Index& mu,
			  typename CubeBasis<IBASIS,DIM>::Support& supp);

  /*!
    For a given wavelet \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <class IBASIS, unsigned int DIM>
  void intersecting_wavelets(const CubeBasis<IBASIS,DIM>& basis,
			     const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename CubeBasis<IBASIS,DIM>::Index>& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
    If this is the case, return true and the intersection of \supp\psi_\lambda
    and \supp\psi_\nu in the form 2^{-j}[k1,k2].
    Otherwise, return false (in this case, j, k1 and k2 will have
    no meaningful values, for performance reasons).
  */
  template <class IBASIS, unsigned int DIM>
  bool intersect_singular_support(const CubeBasis<IBASIS,DIM>& basis,
				  const typename CubeBasis<IBASIS,DIM>::Index& lambda,
				  const typename CubeBasis<IBASIS,DIM>::Index& nu);
}
    
#include <cube/cube_support.cpp>
  
#endif

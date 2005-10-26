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
  /*!
    For a given interval basis IBASIS, compute a cube
      2^{-j}<a,b> = 2^{-j}[a_1,b_1]x...x[a_n,b_n]
    which contains the support of a single primal cube generator
    or wavelet psi_lambda.
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  void support(const CUBEBASIS& basis,
	       const typename CUBEBASIS::Index& lambda,
	       typename CUBEBASIS::Support& supp);

  /*!
    For a given interval basis, compute a cube
     2^{-j}<a,b> = 2^{-j}[a_1,b_1]x...x[a_n,b_n]
    representing the intersection of the support cubes
    corresponding to the indices lambda and mu.
    Function returns true if a nontrivial intersection
    exists, false otherwise. In the latter case 'supp'
    has no meaningful value.
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  bool intersect_supports(const CUBEBASIS& basis,
			  const typename CUBEBASIS::Index& lambda,
			  const typename CUBEBASIS::Index& mu,
			  typename CUBEBASIS::Support& supp);
}
    
#include <cube/cube_support.cpp>
  
#endif

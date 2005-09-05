// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_APPLY_H
#define _WAVELETTL_APPLY_H

#include <algebra/infinite_vector.h>

namespace WaveletTL
{
  /*!
    Apply the stiffness matrix A of an infinite-dimensional equation

      Au = f

    to a given vector v within a specified target accuracy,

      ||w-Av|| <= eta.
      
    The matrix A is assumed to be s^*-compressible, i.e., for a given k,
    there exists a matrix A_k with at most alpha_k*2^k nonzero entries
    per row and column, so that

      ||A-A_k|| <= alpha_k * 2^{-ks}

    for s < s^*.

    References:
    [BDD]  Barinka/Dahlke/Dahmen:
           Adaptive Application of Operators in Standard Representation.
    [CDD1] Cohen/Dahmen/DeVore:
           Adaptive Wavelet Methods for Elliptic Operator Equations - Convergence Rates
  */
  template <class PROBLEM>
  void APPLY(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::WBasis::Index>& v,
	     InfiniteVector<double, typename PROBLEM::WBasis::Index>& w,
	     const double eta);
}

#include <adaptive/apply.cpp>

#endif

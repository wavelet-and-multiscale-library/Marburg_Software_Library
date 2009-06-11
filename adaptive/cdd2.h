// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDD2_H
#define _WAVELETTL_CDD2_H

#include <algebra/infinite_vector.h>

namespace WaveletTL
{
  /*!
    An adaptive solver for the infinite-dimensional problem

      Au = F,

    as developed in [CDD2] and reformulated in [S],[DFR].
    A is assumed to be s.p.d., i.e. the simplifications on [CDD2, p.12] hold.
    Given the problem and a target accuracy \epsilon,
    the algorithm constructs a coefficient vector u_\epsilon, such that
    
      ||u-u_\epsilon|| <= \epsilon.

    The routine has to be given an estimate of ||u|| <= nu = epsilon_0, which may be
    computed beforehand like nu:=||A^{-1}||*||F||.

    References:
    [CDD2] Cohen/Dahmen/DeVore,
           Adaptive Wavelet Methods II - Beyond the Elliptic Case
    [DFR]  Dahlke/Fornasier/Raasch,
           Adaptive Frame Methods for Elliptic Operator Equations
    [S]    Stevenson,
           Adaptive Solution of Operator Equations using Wavelet Frames
  */
  template <class PROBLEM>
  void CDD2_SOLVE(const PROBLEM& P, const double nu, const double epsilon,
		  InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon);
}

#include <adaptive/cdd2.cpp>

#endif

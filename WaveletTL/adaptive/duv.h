// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DUV_H
#define _WAVELETTL_DUV_H

#include <algebra/infinite_vector.h>

namespace WaveletTL
{
  /*
    Some adaptive solvers for the infinite-dimensional problem

      Au = F,

    as developed in [DUV], where A is assumed to be s.p.d.
    Given the problem and a target accuracy \epsilon,
    the algorithm constructs a coefficient vector u_\epsilon, such that
    
      ||u-u_\epsilon|| <= \epsilon.

    The routines has to be given an estimate of ||u|| <= nu = epsilon_0, which may be
    computed beforehand like nu:=||A^{-1}||*||F||.

    References:
    [DUV] Dahmen/Urban/Vorloeper,
          Adaptive Wavelet Methods - Basic Concepts and Applications to the Stokes Problem
  */

  //! steepest descent
  template <class PROBLEM>
  void DUV_SOLVE_SD(const PROBLEM& P, const double nu, const double epsilon,
		    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon, const int jmax = 12);
  
  //! steepest descent quarklet version
  template <class PROBLEM>
  void DUV_QUARKLET_SOLVE_SD(const PROBLEM& P, const double nu, const double epsilon,
		    InfiniteVector<double, typename PROBLEM::QuarkletFrame::Index>& u_epsilon, CompressionStrategy strategy, const int pmax = 0, 
                    const int jmax = 12, const double a = 2, const double b = 2);
}

#include <adaptive/duv.cpp>

#endif

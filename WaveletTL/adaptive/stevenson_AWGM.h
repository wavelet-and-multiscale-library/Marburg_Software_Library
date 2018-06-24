// -*- c++ -*-

// +--------------------------------------------------------------------+
// | stevenson_AWGM.h, Copyright (c) 2018                               |
// | Henning Zickermann <zickermann@mathematik.uni-marburg.de>          |
// |                                                                    |
// | This file is part of WaveletTL - the Wavelet Template Library.     |
// |                                                                    |
// | Contact: AG Numerik, Philipps University Marburg                   |
// |          http://www.mathematik.uni-marburg.de/~numerik/            |
// +--------------------------------------------------------------------+


#ifndef _WAVELETTL_STEVENSON_AWGM_H
#define _WAVELETTL_STEVENSON_AWGM_H

#include <set>
#include <algebra/infinite_vector.h>
#include <adaptive/compression.h>
#include <adaptive/apply.h>
#include <utils/convergence_logger.h>


namespace WaveletTL
{


/*
  An optimal, adaptive Wavelet-Galerkin method (AWGM) without coarsening of the iterands
  as developed in [GHS07].
  The algorithm applies to linear operator equations, reformulated as infinite-dimensional
  matrix-vector equation

   Au = F

  in \ell_2 by means of a Wavelet basis, where A is assumed to be boundedly invertible,
  symmetric and positive-definite.
  Given the problem and a target accuracy epsilon, the algorithm constructs a coefficient vector
  u_epsilon, such that the \ell_2-norm of the residual is lesser than or equal to epsilon, i.e.

    ||F-Au_epsilon||_2 <= epsilon.

  You can specify a maximal level jmax for the internal APPLY calls.

  References:
  [GHS07]  T. Gantumur, H. Harbrecht, R.P. Stevenson, An Optimal Adaptive Wavelet Method
           without Coarsening of the Iterands, Math. Comp., 76:615–629, 2007.

  [Ste09]  R.P. Stevenson, Adaptive wavelet methods for solving operator equations:
           An overview, Multiscale, Nonlinear and Adaptive Approximation: 543-597.
           Springer-Verlag Berlin Heidelberg, 2009.
*/



using std::set;
using MathTL::InfiniteVector;



/*
 * The routine SOLVE from [GHS07] with parameters alpha, omega, gamma, theta > 0.
 * In [GHS07], SOLVE was proven to be of optimal computational complexity in case 0 < omega < alpha < 1,
 * (alpha + omega)/(1-omega) < kappa(A)^{-1/2} and 0 < gamma < 1/6* kappa(A)^{-1/2}*(alpha-omega)/(1+omega).
 * However, in practice a better performance can be reached when choosing the parameters outside these ranges.
 */
template <class PROBLEM>
void AWGM_SOLVE(const PROBLEM& P, const double epsilon,
                InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
                const int jmax,
                MathTL::AbstractConvergenceLogger& logger = MathTL::DummyLogger(),
                const double alpha = 0.9,
                const double omega = 0.01,
                const double gamma = 0.01,
                const double theta = 0.4,
                const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess = InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>(),
  #if _WAVELETTL_USE_TBASIS == 1
                const CompressionStrategy strategy = tensor_simple);
  #else
                const CompressionStrategy strategy = St04a);
  #endif



/*
 * The routine SOLVE from [GHS07] with additional possibility to specify nu_{-1}.
 */
template <class PROBLEM>
void AWGM_SOLVE(const PROBLEM& P, const double epsilon,
                InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
                const int jmax,
                const double nu_neg1,
                MathTL::AbstractConvergenceLogger& logger = MathTL::DummyLogger(),
                const double alpha = 0.9,
                const double omega = 0.01,
                const double gamma = 0.01,
                const double theta = 0.4,
                const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess = InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>(),
  #if _WAVELETTL_USE_TBASIS == 1
                const CompressionStrategy strategy = tensor_simple);
  #else
                const CompressionStrategy strategy = St04a);
  #endif



/*
 * A simplified version of GALSOLVE from [GHS07].
 */
template <class PROBLEM>
void GALSOLVE(const PROBLEM& P, const set<typename PROBLEM::WaveletBasis::Index>& Lambda,
              const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& g_Lambda,
              InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& w_Lambda,
              const double delta,
              const double epsilon);

}


#include "stevenson_AWGM.cpp"

#endif // _WAVELETTL_STEVENSON_AWGM_H

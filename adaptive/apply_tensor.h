// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_APPLY_TENSOR_H
#define	_WAVELETTL_APPLY_TENSOR_H

#include <algebra/infinite_vector.h>
#include <adaptive/compression.h>
#include <utils/array1d.h>
#include <utils/tiny_tools.h>

#include <iostream>
#include <list>
#include <map>
#include <limits> // to use -inf as a literal
#include <algorithm>

using MathTL::Array1D;


namespace WaveletTL
{
    /*!
    These APPLY routines are tailored for the case that the
    underlying wavelet basis has a vector valued minimal level j, e.g.
    bases of tensor structure type as modeled by tbasis/qtbasis.

    Apply the stiffness matrix A of an infinite-dimensional equation

      Au = f

    to a given vector v within a specified target accuracy,

      ||w-Av|| <= eta.

    The matrix A is assumed to be s^*-compressible, i.e., for a given k,
    there exists a matrix A_k with at most C*k nonzero entries
    per row and column, so that

      ||A-A_k|| <= alpha_k * k^{-s} <= alpha * k^{-s}  if s^*< \infty,
      ||A-A_k|| <= alpha_k * 2^{-rho*k} <= alpha * 2^{-rho*k} if s^* = \infty.

    The basis provided by tbasis gives s^*=\infty for operator equations with 
    bounded coefficients. And uses C*k^dim many entries per column in A_k.
    Right now only the case s^*=\infty is implemented.

    The class PROBLEM has to provide the routine

      double alphak(const unsigned int) const

    that is used as the parameter alpha (independent of k)

    An example of the template parameter PROBLEM is the template class
    SturmEquation<WBASIS>.

    The action of APPLY can be restricted to entries below a maximal scale jmax.

    References:
    [DSS]  Dijkema, Schwab, Stevenson:
           An Adaptive Wavelet Method for Solving High-Dimensional Elliptic PDEs
  */

  /*
   * Implmentation of APPLY from [DSS]
   * 
   * The coeffs of the vector v are sorted into buckets according to their absolute value. 
   * A different approximation of A is applied to each bucket.
   * The optimal approximation matrix for each bucket is computed directly. 
   * The formula for L2-orthogonal wavelets from [DSS] is used. 
   * This may lead to an overestimate and thus too high computational effort.
   * 
   * In a previous version of the code the optimal effort was computed by a costly minimization problem (-> commented code)
   */
  template <class PROBLEM>
  void APPLY_TENSOR(PROBLEM& P,
          const InfiniteVector<double, typename PROBLEM::Index>& v,
          const double eta,
          InfiniteVector<double, typename PROBLEM::Index>& w,
          const int jmax = 99,
          const CompressionStrategy strategy = tensor_simple,
          const bool preconditioning = true);
  
  template <class PROBLEM>
  void APPLY_TENSOR(PROBLEM& P,
          const InfiniteVector<double, int>& v,
          const double eta,
          InfiniteVector<double, int>& w,
          const int jmax = 99,
          const CompressionStrategy strategy = tensor_simple,
          const bool preconditioning = true);

#if 0
  /* discontinued code:
   * 
   * In order to compute the approximation of Av, v is sorted into buckets (v_[p]).
   * The entries in each bucket are applied to different approximations A^Jp of A.
   * Better approximations A^jp are more costly. Therefore the optimal compression/benefit ratio
   * has to be computed. cost((j_p))-> min; approximation_error <= eta
   * The minimiation problem may use an inital guess for (j_p)
   * 
   * There is commented code for a restriction to levels given by a levelvector
   */
  template <class PROBLEM>
  void APPLY(const PROBLEM& P,
          const InfiniteVector<double, typename PROBLEM::Index>& v,
          const double eta,
          //const set<unsigned int> levelwindow,
          Array1D<int>& jp_guess,
          InfiniteVector<double, typename PROBLEM::Index>& w,
          const int jmax = 99,
          const CompressionStrategy strategy = tensor_simple,
          const bool preconditioning = true);
#endif
}

#include <adaptive/apply_tensor.cpp>

#endif	/* _WAVELETTL_APPLY_TENSOR_H */


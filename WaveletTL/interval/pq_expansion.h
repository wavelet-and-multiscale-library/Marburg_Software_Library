// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PQ_EXPANSION_H
#define _WAVELETTL_PQ_EXPANSION_H

#include <algebra/infinite_vector.h>
#include <interval/pq_frame.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT> class PQFrame;

  /*!
    helper function, integrate a smooth function f against a
    primal [P] generator or wavelet
  */
  template <int d, int dT>
  double integrate(const Function<1>* f,
		   const PQFrame<d,dT>& basis,
		   const typename PQFrame<d,dT>::Index& lambda);

  /*!
    helper function, integrate two primal [P] generators or wavelets
    against each other (for the Gramian)
  */
  template <int d, int dT>
  double integrate(const PQFrame<d,dT>& basis,
		   const typename PQFrame<d,dT>::Index& lambda,
		   const typename PQFrame<d,dT>::Index& mu);
  
  /*!
    For a given function, compute all integrals w.r.t. the primal
    or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
    - When integrating against the primal functions, the integrand has to be smooth
      to be accurately reproduced by the dual basis.
    - When integration against dual functions is specified,
      we integrate against the primal ones instead and multiply the resulting
      coefficients with the inverse of the primal gramian.

    Maybe a thresholding of the returned coefficients is helpful (e.g. for
    expansions of spline functions).
  */
  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const PQFrame<d,dT>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename PQFrame<d,dT>::Index>& coeffs,
         const int pmax = 0);

  /*!
    analogous routine for Vector<double> output
  */
  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const PQFrame<d,dT>& basis,
	 const bool primal,
	 const int jmax,
	 Vector<double>& coeffs);
  
  template <int d, int dT>
  double factor(const PQFrame<d,dT>& basis, const int l, const int k, const int r, const int q, const int j);
  
  template <int d, int dT>
  void setup_factor_matrix(const PQFrame<d,dT>& basis, SparseMatrix<double>& A_Lambda, const int j, const int p);
}

#include <interval/pq_expansion.cpp>

#endif

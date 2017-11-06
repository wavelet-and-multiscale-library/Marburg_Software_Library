// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_GALERKIN_UTILS_H
#define _WAVELETTL_GALERKIN_UTILS_H

#include <set>

#include <algebra/sparse_matrix.h>
#include <algebra/vector.h>
#include <omp.h>

using MathTL::SparseMatrix;
using MathTL::Vector;

namespace WaveletTL
{
  /*!
    Setup the (sparse, preconditioned per default) stiffness matrix for a given problem and a given active
    index set Lambda.
   * Integer version should be faster
  */
  template <class PROBLEM>
  void setup_stiffness_matrix(PROBLEM& P,
			      const std::set<typename PROBLEM::Index>& Lambda,
			      SparseMatrix<double>& A_Lambda,
			      bool preconditioned = true); 

    /*!
     * Same, but with integers instead of wavelet indices (faster)
     * PROBLEM& P is not const because then you are able to use also cached problems
     */
    template <class PROBLEM>
    void setup_stiffness_matrix(PROBLEM& P,
            const std::set<int>& Lambda,
            SparseMatrix<double>& A_Lambda,
            bool preconditioned = true);
  /*!
    Setup the (sparse, preconditioned per default) stiffness matrix for a given problem and a given active
    index set Lambda.
  */
  template <class PROBLEM>
  void setup_stiffness_matrix(PROBLEM& P,
			      const std::set<typename PROBLEM::Index>& Lambda1,
			      const std::set<typename PROBLEM::Index>& Lambda2,
			      SparseMatrix<double>& A_Lambda,
			      bool preconditioned = true); 


  
  /*!
   * Setup the (preconditioned) right-hand side for a given problem and a given active
   * index set Lambda.
   * 
   * CAUTION: Index and int version of this method yield different results!
   * 
   * Index version of the method assumes that P.f and P.D exist and returns (P.f(lambda)/P.D(lambda))_lambda
   * Integer Version assumes that P.f already returns the preconditioned coeffs, i.e, there is no call to P.D
  */
  template <class PROBLEM>
  void setup_righthand_side(PROBLEM& P,
			    const std::set<typename PROBLEM::Index>& Lambda,
			    Vector<double>& F_Lambda);
    template <class PROBLEM>
    void setup_righthand_side(PROBLEM& P,
            const std::set<int>& Lambda,
            Vector<double>& F_Lambda);
}

#include <galerkin/galerkin_utils.cpp>

#endif

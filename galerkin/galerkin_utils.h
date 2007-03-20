// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_GALERKIN_UTILS_H
#define _WAVELETTL_GALERKIN_UTILS_H

#include <set>

#include <algebra/sparse_matrix.h>
#include <algebra/vector.h>

using MathTL::SparseMatrix;
using MathTL::Vector;

namespace WaveletTL
{
  /*!
    Setup the (sparse, preconditioned per default) stiffness matrix for a given problem and a given active
    index set Lambda.
  */
  template <class PROBLEM>
  void setup_stiffness_matrix(const PROBLEM& P,
			      const std::set<typename PROBLEM::Index>& Lambda,
			      SparseMatrix<double>& A_Lambda,
			      bool preconditioned = true); 
  /*!
    Setup the (sparse, preconditioned per default) stiffness matrix for a given problem and a given active
    index set Lambda.
  */
  template <class PROBLEM>
  void setup_stiffness_matrix(const PROBLEM& P,
			      const std::set<typename PROBLEM::Index>& Lambda1,
			      const std::set<typename PROBLEM::Index>& Lambda2,
			      SparseMatrix<double>& A_Lambda,
			      bool preconditioned = true); 


  
  /*!
    Setup the (preconditioned) right-hand side for a given problem and a given active
    index set Lambda.
  */
  template <class PROBLEM>
  void setup_righthand_side(const PROBLEM& P,
			    const std::set<typename PROBLEM::Index>& Lambda,
			    Vector<double>& F_Lambda); 
}

#include <galerkin/galerkin_utils.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PRECOND_H
#define _WAVELETTL_PRECOND_H

#include <algebra/infinite_vector.h>
#include <algebra/infinite_matrix.h>

using MathTL::InfiniteVector;
using MathTL::InfiniteDiagonalMatrix;

namespace WaveletTL
{
  /*!
    Base class for preconditioners in wavelet-Galerkin schemes
   */
  template <class PROBLEM>
  class WaveletGalerkinPreconditioner
  {
  public:
    /*!
      apply the preconditioner to a given vector y, i.e., solve Px = y
    */
    virtual void apply_preconditioner(const InfiniteVector<double, typename PROBLEM::Index>& y,
				      InfiniteVector<double, typename PROBLEM::Index>& x) const = 0;
  };

  /*!
    Base class for fully diagonal preconditioners P=D^{-1}, D=diag(d_i)
  */
  template <class PROBLEM>
  class DiagonalPreconditioner
    : public WaveletGalerkinPreconditioner<PROBLEM>,
      public InfiniteDiagonalMatrix<double, typename PROBLEM::Index>
  {
  public:
    /*!
      evaluate the diagonal preconditioner D
    */
    virtual double diag(const typename PROBLEM::Index& lambda) const = 0;

    /*!
      apply preconditioner
    */
    void apply_preconditioner(const InfiniteVector<double, typename PROBLEM::Index>& y,
			      InfiniteVector<double, typename PROBLEM::Index>& x) const;
  };
}

#include <galerkin/precond.cpp>

#endif

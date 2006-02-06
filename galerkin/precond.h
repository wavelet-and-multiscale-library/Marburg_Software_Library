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
      apply preconditioner, x = D^{-1}y
    */
    void apply_preconditioner(const InfiniteVector<double, typename PROBLEM::Index>& y,
			      InfiniteVector<double, typename PROBLEM::Index>& x) const;
  };

  /*!
    Fully diagonal preconditioner based on (NE) for the energy space, i.e.,
      D=diag(2^{t|lambda|})
  */
  template <class PROBLEM>
  class DiagonalDyadicPreconditioner
    : public DiagonalPreconditioner<PROBLEM>
  {
  public:
    /*!
      default constructor, takes the exponent t
    */
    DiagonalDyadicPreconditioner(const double expo) : t(expo) {}

    /*!
      evaluate the diagonal preconditioner D
    */
    double diag(const typename PROBLEM::Index& lambda) const { return ldexp(1.0, t*lambda.j()); }

  protected:
    //! exponent t
    double t;
  };

//   /*!
//     Fully diagonal preconditioner using the energy norms D=sqrt(diag(A)), i.e.,
//     D depends on a given unpreconditioned problem
//   */
//   template <class PROBLEM>
//   class DiagonalEnergyNormPreconditioner
//     : public DiagonalPreconditioner<PROBLEM>
//   {
//   public:
//     /*!
//       default constructor, takes a given problem
//     */
//     DiagonalEnergyNormPreconditioner(const PROBLEM* p) : problem(p) {}
    
//     /*!
//       evaluate the diagonal preconditioner D
//     */
//     double diag(const typename PROBLEM::Index& lambda) const { return sqrt(problem.a(lambda,lambda)); }
    
//   protected:
//     //! an instance of the problem
//     const PROBLEM* problem;
//   };
}

#include <galerkin/precond.cpp>

#endif

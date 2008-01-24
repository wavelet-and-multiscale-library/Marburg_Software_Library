// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_HELMHOLTZ_EQUATION_H
#define _WAVELETTL_HELMHOLTZ_EQUATION_H

#include <utils/function.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/full_helmholtz.h>
#include <galerkin/infinite_preconditioner.h>
#include <galerkin/gramian.h>
#include <galerkin/poisson_equation.h>
#include <galerkin/cached_problem.h>
#include <adaptive/compression.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    This class models the wavelet discretization of the one--dimensional
    Helmholtz equation

      -u''(t) + alpha * u(t) = f(t), 0 <= t <= 1

    with Dirichlet boundary conditions u(0)=u(1)=0.

    The wavelet bases useable for the discretization are those
    modeled by SplineBasis<d,dT>.

    Internally, HelmholtzEquation1D holds two distinct caches for the
    Gramian (identity) part of the stiffness matrix and the Laplacian part.
  */
  template <int d, int dT, int J0>
  class HelmholtzEquation1D
    :  public FullyDiagonalEnergyNormPreconditioner<typename SplineBasis<d,dT,P_construction,1,1,0,0,J0>::Index>
  {
  public:
    /*!
      type of the wavelet basis
    */
    typedef SplineBasis<d,dT,P_construction,1,1,0,0,J0> WaveletBasis;

    /*!
      constructor from a given wavelet basis and a right-hand side y
    */
    HelmholtzEquation1D(const WaveletBasis& basis,
			const double alpha,
			const InfiniteVector<double, typename WaveletBasis::Index>& y);

    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return basis_; }
    
    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

    /*!
      index type of vectors and matrices
    */
    typedef typename Vector<double>::size_type size_type;

    /*!
      space dimension of the problem
    */
    static const int space_dimension = 1;

    /*!
      differential operators are local
    */
    static bool local_operator() { return true; }

    /*!
      (half) order t of the operator
      (inherited from FullyDiagonalDyadicPreconditioner)
    */
    double operator_order() const { return 1.; }
    
    /*!
      set (unpreconditioned) right-hand side y,
      call this _after_ the corresponding set_alpha()
    */
    void set_rhs(const InfiniteVector<double,Index>& y) const;
    
    /*!
      set reaction coefficient alpha
    */
    void set_alpha(const double alpha) const;

    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const Index& lambda) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a
      (inherited from FullyDiagonalEnergyNormPreconditioner)
     */
    double a(const Index& lambda,
 	     const Index& nu) const;

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta, InfiniteVector<double,Index>& coeffs) const
    {
      coeffs = y_precond; // dirty
    }

    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const {
      return y_.get_coefficient(lambda);
    }

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return l2_norm(y_); }

    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const;

    /*!
      estimate the spectral norm ||A^{-1}||
    */
    double norm_Ainv() const;

    /*!
      estimate compressibility exponent s^*
    */
    double s_star() const {
      return 1.0 + WaveletBasis::primal_vanishing_moments(); // [St04a], Th. 2.3 for n=1
    }

    /*!
      estimate the compression constants alpha_k in
        ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 2*norm_A(); // suboptimal
    }

    /*!
      w += factor * (stiffness matrix entries in column lambda on level j)
    */
    void add_level (const Index& lambda,
		    //InfiniteVector<double, Index>& w,
		    Vector<double>& w,
		    const int j,
		    const double factor,
		    const int J,
		    const CompressionStrategy strategy = St04a) const;

  protected:
    WaveletBasis basis_;
    mutable double alpha_;
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_;
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_precond;

    mutable FullHelmholtz<d,dT,J0> H_;

  public:
    IntervalGramian<WaveletBasis> G_;
    CachedProblem<IntervalGramian<WaveletBasis> > GC_;
    PoissonEquation1D<d,dT,J0> A_;
    CachedProblem<PoissonEquation1D<d,dT,J0> > AC_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };

}

#include <galerkin/helmholtz_equation.cpp>

#endif

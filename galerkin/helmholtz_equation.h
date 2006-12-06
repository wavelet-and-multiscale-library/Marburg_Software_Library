// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
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
  */
  template <int d, int dT>
  class HelmholtzEquation1D
    :  public FullyDiagonalEnergyNormPreconditioner<typename SplineBasis<d,dT>::Index>
  {
  public:
    /*!
      constructor from a right-hand side f
    */
    HelmholtzEquation1D(const Function<1>* f,
			const double alpha = 0,
			const bool precompute_rhs = true);

    /*!
      type of the wavelet basis
    */
    typedef SplineBasis<d,dT> WaveletBasis;

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
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite Gauss quadrature rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
    double a(const Index& lambda,
 	     const Index& nu,
 	     const unsigned int p) const;

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta, InfiniteVector<double,Index>& coeffs) const;

    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const;

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }

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
		    InfiniteVector<double, Index>& w, const int j,
		    const double factor,
		    const int J,
		    const CompressionStrategy strategy = St04a) const;

  protected:
    const Function<1>* f_;
    double alpha_;
    SplineBasis<d,dT> basis_;
    FullHelmholtz<d,dT> A_;

    // flag whether right-hand side has already been precomputed
    mutable bool rhs_precomputed;

    /*!
      precomputation of the right-hand side
      (constness is not nice but necessary to have RHS a const function)
    */
    void precompute_rhs() const;

    // (preconditioned) right-hand side coefficients on a fine level
    mutable InfiniteVector<double,Index> fcoeffs_unsorted;
    
    // right-hand side coefficients on a fine level, sorted by modulus
    mutable Array1D<std::pair<Index,double> > fcoeffs;

    // (squared) \ell_2 norm of the precomputed right-hand side
    mutable double fnorm_sqr;

    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };



}

#include <galerkin/helmholtz_equation.cpp>

#endif

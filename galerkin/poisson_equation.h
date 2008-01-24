// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_POISSON_EQUATION_H
#define _WAVELETTL_POISSON_EQUATION_H

#include <utils/function.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/full_laplacian.h>
#include <galerkin/infinite_preconditioner.h>
#include <adaptive/compression.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    The Poisson equation -u''=f, u(0)=u(1)=0 in one space dimension.
  */
  template <int d, int dT, int J0>
  class PoissonEquation1D
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
    PoissonEquation1D(const WaveletBasis& basis,
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
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const typename WaveletBasis::Index& lambda) const {
      return y_.get_coefficient(lambda);
    }

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
 	     InfiniteVector<double,typename WaveletBasis::Index>& coeffs) const {
      coeffs = y_; // dirty
      coeffs.scale(this, -1); // coeffs *= D^{-1}
    }

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return l2_norm(y_); }

    /*!
      set right-hand side y
    */
    void set_rhs(const InfiniteVector<double,Index>& y) const {
      y_ = y;
    }

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
    const WaveletBasis& basis_;
    mutable FullLaplacian<d,dT,J0> A_;
    
    // rhs, mutable to have 'const' method
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };

}

#include <galerkin/poisson_equation.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_RING_HELMHOLTZ_H
#define _WAVELETTL_RING_HELMHOLTZ_H

#include <utils/function.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <interval/spline_basis.h>
#include <ring/ring_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <galerkin/ring_gramian.h>
#include <galerkin/ring_laplacian.h>
#include <galerkin/cached_problem.h>
// #include <adaptive/compression.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    This class models the Helmholtz equation -Delta u+alpha*u=f
    over a ring-shaped domain, to be used in adaptive algorithms.

    Internally, RingHelmholtzEquation holds two distinct caches for the
    Gramian (identity) part of the stiffness matrix and the Laplacian part.
  */
  template <int d, int dt, int s0, int s1>
  class RingHelmholtzEquation
    : public FullyDiagonalEnergyNormPreconditioner<typename RingBasis<d,dt,s0,s1>::Index>
  {
  public:
    //! wavelet basis class
    typedef RingBasis<d,dt,s0,s1> WaveletBasis;

    //! wavelet index class
    typedef typename WaveletBasis::Index Index;
    
    //! type of 1D basis in angular direction
    typedef typename RingBasis<d,dt,s0,s1>::Basis0 Basis0;
    
    //! type of 1D indices in angular direction
    typedef typename Basis0::Index Index0;
    
    //! type of 1D basis in radial direction
    typedef typename RingBasis<d,dt,s0,s1>::Basis1 Basis1;

    //! type of 1D indices in radial direction
    typedef typename Basis1::Index Index1;
  
    //! index type of vectors and matrices
    typedef typename Vector<double>::size_type size_type;

    /*!
      constructor from a given wavelet basis,
      two filenames for the precomputed Gramian and Laplacian (plus corresponding jmax)
      and a right-hand side y
    */
    RingHelmholtzEquation(const WaveletBasis& basis,
			  const char* G_file,
			  const char* A_file,
			  const int jmax,
			  const double alpha,
			  const InfiniteVector<double, typename WaveletBasis::Index>& y);

    //! read access to the basis
    const WaveletBasis& basis() const { return basis_; }
    
    //! space dimension of the problem
    static const int space_dimension = 2;

    //! identity operator is local
    static bool local_operator() { return true; }

    //! (half) order t of the operator
    static double operator_order() { return 1; }

    //! evaluate the diagonal preconditioner D
    double D(const Index& lambda) const { return sqrt(a(lambda,lambda)); }
    
    //! evaluate the (unpreconditioned) bilinear form a
    double a(const Index& lambda,
 	     const Index& mu) const;

    //! estimate the spectral norm ||A||
    double norm_A() const;

    //! estimate the spectral norm ||A^{-1}||
    double norm_Ainv() const;

    //! estimate compressibility exponent s^*
    double s_star() const {
      return WaveletBasis::primal_vanishing_moments();
    }
    
    /*!
      estimate the compression constants alpha_k in
      ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 2*norm_A(); // suboptimal
    }

    //! evaluate the (unpreconditioned) right-hand side f
    double f(const Index& lambda) const {
      return y_.get_coefficient(lambda);
    }

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
  	     InfiniteVector<double,Index>& coeffs) const {
      coeffs = y_; // dirty
    }

    //! compute (or estimate) ||F||_2
    double F_norm() const { return l2_norm(y_); }

    //! set right-hand side y
    void set_rhs(const InfiniteVector<double,Index>& y) const {
      y_ = y;
    }

    /*!
      set reaction coefficient alpha
    */
    void set_alpha(const double alpha) const;

  protected:
    const WaveletBasis basis_;
    mutable double alpha_;
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_;
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_precond;

  public:
    RingGramian<d,dt,s0,s1> G_;
    CachedProblemFromFile<RingGramian<d,dt,s0,s1> > GC_;
    RingLaplacian<d,dt,s0,s1> A_;
    CachedProblemFromFile<RingLaplacian<d,dt,s0,s1> > AC_;

    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };

}

#include <galerkin/ring_helmholtz.cpp>

#endif

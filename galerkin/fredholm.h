// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_FREDHOLM_H
#define _WAVELETTL_FREDHOLM_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <algebra/vector.h>
#include <interval/spline_basis.h>
#include <galerkin/cached_problem.h>
#include <adaptive/compression.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    This base class models the system matrix K^*K of a zero-order Fredholm integral operator

      Au(t) = int_0^1 k(s,t)u(s) ds

    with respect to a given wavelet basis on the interval [0,1], i.e.,
 
      <A^*Au,v> = int_0^1 (int_0^1 k(s,t)u(s) ds)(int_0^1 k(s,t)v(s) ds) dt
                = int_0^1 int_0^1 g(s,s) u(s) v(t) ds dt,

      g(s,t)    = int_0^1 k(s,x) k(t,x) dt.

    The kernel g should be known explicitily. This is the case, e.g., for the Volterra
    integration operator

      k(s,t) = / 1 for 0<=s<=t
             = \ 0 otherwise

    with

      g(s,t) = 1-max(s,t).

    The class has the minimal signature to be used within the APPLY routine
    or within adaptive solvers like CDD1.
  */
  template <int d, int dT, int J0>
  class FredholmIntegralOperator
    :  public FullyDiagonalEnergyNormPreconditioner<typename SplineBasis<d,dT,P_construction,0,0,0,0,J0>::Index>
  {
  public:
    /*!
      type of the wavelet basis
    */
    typedef SplineBasis<d,dT,P_construction,0,0,0,0,J0> WaveletBasis;
    
    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;
    
    /*!
      index type of vectors and matrices
    */
    typedef typename Vector<double>::size_type size_type;

    /*!
      constructor from a given wavelet basis and a given right-hand side y
    */
    FredholmIntegralOperator(const WaveletBasis& basis,
			     const InfiniteVector<double,Index>& y);

    /*!
      purely virtual destructor
    */
    virtual ~FredholmIntegralOperator() = 0;
    
    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return basis_; }

    /*!
      space dimension of the problem
    */
    static const int space_dimension = 1;
    
    /*!
      Fredholm operators are in general nonlocal
    */
    static bool local_operator() { return false; }
    
    /*!
      (half) order t of the operator
    */
    static double operator_order() { return 0; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const Index& lambda) const { return 1; }

    /*
      kernel function
    */
    virtual double g(const double s, const double t) const = 0;
    
    /*!
      evaluate the (unpreconditioned) bilinear form a
    */
    double a(const Index& lambda,
 	     const Index& nu) const;
    
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
      return WaveletBasis::primal_vanishing_moments();
    }
     
    /*!
      estimate the compression constants alpha_k in
      ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 2*norm_A(); // suboptimal
    }
    
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
    }
    
    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return l2_norm(y_); }
    
    /*!
      set right-hand side y
    */
    void set_rhs(const InfiniteVector<double, typename WaveletBasis::Index>& y) const {
      y_ = y;
    }
 
  protected:
    // the wavelet basis
    const WaveletBasis& basis_;
    
    // rhs, mutable to have 'const' method
    mutable InfiniteVector<double,Index> y_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };

  /*!
    example: Volterra integration operator
  */
  template <int d, int dT, int J0>
  class VolterraIntegralOperator
    :  public FredholmIntegralOperator<d,dT,J0>
  {
  public:
    /*!
      type of the wavelet basis
    */
    typedef SplineBasis<d,dT,P_construction,0,0,0,0,J0> WaveletBasis;
    
    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;
    
    /*!
      index type of vectors and matrices
    */
    typedef typename Vector<double>::size_type size_type;

    /*!
      constructor from a given wavelet basis and a given right-hand side y
    */
    VolterraIntegralOperator(const WaveletBasis& basis,
			     const InfiniteVector<double,Index>& y);
    
    /*
      kernel function
    */
    double g(const double s, const double t) const {
      return 1.0 - std::max(s,t);
    }
  };
}

#include <galerkin/fredholm.cpp>

#endif

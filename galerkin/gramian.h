// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_GRAMIAN_H
#define _WAVELETTL_GRAMIAN_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <algebra/vector.h>
#include <interval/spline_basis.h>
#include <galerkin/full_gramian.h>
#include <adaptive/compression.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    This class models equations Ax=y with the biinfinite Gramian matrix 

      A = <Psi,Psi>

    of a given (DS-/DKU-) wavelet basis on the unit interval [0,1],
    to be used in (adaptive) algorithms.

    The class has the minimal signature to be used within the APPLY routine
    or within adaptive solvers like CDD1.
  */
  template <class WBASIS>
  class IntervalGramian
  {
  public:
    /*!
      constructor from a given wavelet basis and a given right-hand side y
    */
    IntervalGramian(const WBASIS& basis,
	    const InfiniteVector<double, typename WBASIS::Index>& y);
    
    /*!
      make template argument accessible
    */
    typedef WBASIS WaveletBasis;

    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return basis_; }
    
    /*!
      index type of vectors and matrices
    */
    typedef typename Vector<double>::size_type size_type;

    /*!
      space dimension of the problem
    */
    static const int space_dimension = 1;

    /*!
      identity operator is local
    */
    static bool local_operator() { return true; }

    /*!
      (half) order t of the operator
    */
    static double operator_order() { return 0; }

    /*!
      evaluate the diagonal preconditioner D (essentially, we don't need any)
    */
    double D(const typename WaveletBasis::Index& lambda) const { return 1; }

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
    const WaveletBasis& basis_;
    
    // rhs, mutable to have 'const' method
    mutable InfiniteVector<double, typename WBASIS::Index> y_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };

  
  /*!
    template specialization to SplineBasis
  */
  template <int d, int dT, int s0, int s1, int sT0, int sT1>
  class IntervalGramian<SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1> >
  {
  public:
    /*!
      make template argument accessible
    */
    typedef SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1> WaveletBasis;
    
    /*!
      constructor from a given wavelet basis and a given right-hand side y
    */
    IntervalGramian(const SplineBasis<d,dT,P_construction,s0,s1,sT0,sT1>& basis,
 		    const InfiniteVector<double, typename WaveletBasis::Index>& y);
    
    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;
    
    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return basis_; }
    
    /*!
      index type of vectors and matrices
    */
    typedef typename Vector<double>::size_type size_type;

    /*!
      space dimension of the problem
    */
    static const int space_dimension = 1;

    /*!
      identity operator is local
    */
    static bool local_operator() { return true; }

    /*!
      (half) order t of the operator
    */
    static double operator_order() { return 0; }

    /*!
      evaluate the diagonal preconditioner D (essentially, we don't need any)
    */
    double D(const typename WaveletBasis::Index& lambda) const { return 1; }

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
    mutable FullGramian<d,dT,s0,s1,sT0,sT1> G_;
    
    // rhs, mutable to have 'const' method
    mutable InfiniteVector<double, typename WaveletBasis::Index> y_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };

}

#include <galerkin/gramian.cpp>

#endif

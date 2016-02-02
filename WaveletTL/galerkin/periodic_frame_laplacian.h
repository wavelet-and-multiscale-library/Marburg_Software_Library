// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PERIODIC_FRAME_LAPLACIAN_H
#define _WAVELETTL_PERIODIC_FRAME_LAPLACIAN_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <algebra/vector.h>
#include <interval/periodic.h>
#include <numerics/periodiclp.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    This class models equations Ax=y with the biinfinite Laplacian matrix 

      A = <Psi,Psi>

    of a given periodic wavelet basis on the unit interval [0,1],
    to be used in (adaptive) algorithms.

    The class has the minimal signature to be used within the APPLY routine
    or within adaptive solvers like CDD1.
  */
  template <class RFRAME>
  class PeriodicFrameIntervalLaplacian
  {
  public:
    /*!
      constructor from a given wavelet basis and a given right-hand side y
    */
    PeriodicFrameIntervalLaplacian(const PeriodicLaplacianProblem& plp, const PeriodicFrame<RFRAME>& frame);
    
    /*!
      make template argument accessible
    */
    typedef PeriodicFrame<RFRAME> WaveletBasis;

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
    static double operator_order() 
        { 
            return 1;
        }

    /*!
      evaluate the diagonal preconditioner D (dyadic scaling)
    */
    double D(const typename WaveletBasis::Index& lambda) const 
    { 
        
        return ldexp(1.0, lambda.j())*pow(1+lambda.p(),5); //2^j*(p+1)^2+delta (delta=3)
        //return ldexp(1.0, lambda.j()); //2^j
        //return 1;
        //return sqrt(a(lambda,lambda));
    }

    /*!
      evaluate the (unpreconditioned) bilinear form a
    */
    double a(const typename WaveletBasis::Index& lambda,
  	     const typename WaveletBasis::Index& nu) const;

    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const {
      return normA;
    }

    /*!
      estimate the spectral norm ||A^{-1}||
    */
    double norm_Ainv() const {
      return normAinv;
    }

    /*!
      estimate compressibility exponent s^*
    */
    double s_star() const {
      return WaveletBasis::primal_vanishing_moments();
      //return 0.7;
      //return 0.5;
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
    /*double f(const typename WaveletBasis::Index& lambda) const {
      return y_.get_coefficient(lambda);
    }*/
    double f(const Index& lambda) const;

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
 	     InfiniteVector<double, Index>& coeffs) const; /*{
      coeffs = y_; // dirty
    }*/

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }

    /*!
      set right-hand side y
    */
    /*void set_rhs(const InfiniteVector<double, typename WaveletBasis::Index>& y) const {
      y_ = y;
    }*/

  protected:
    const PeriodicLaplacianProblem& plp_;
    WaveletBasis basis_;
    
    
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
    
    // flag whether right-hand side has already been precomputed
    mutable bool rhs_precomputed;
    
    // right-hand side coefficients on a fine level, sorted by modulus
    mutable Array1D<std::pair<Index,double> > fcoeffs;
    
    // (squared) \ell_2 norm of the precomputed right-hand side
    mutable double fnorm_sqr;
    
    /*!
      precomputation of the right-hand side
      (constness is not nice but necessary to have RHS a const function)
    */
    void precompute_rhs() const;
  };

  
#if 0
  class PeriodicLaplacianProblem
  {
  public:
    /*!
      virtual destructor
     */
    virtual ~PeriodicLaplacianProblem ();
    virtual double g(const double t) const = 0;
  };
#endif
}

#include <galerkin/periodic_frame_laplacian.cpp>

#endif

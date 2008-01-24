// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_RING_LAPLACIAN_H
#define _WAVELETTL_RING_LAPLACIAN_H

#include <set>
#include <map>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>
#include <algebra/vector.h>
#include <adaptive/compression.h>
#include <ring/ring_basis.h>
#include <galerkin/infinite_preconditioner.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    This class models the "polar coordinate" Laplacian matrix
    over a ring-shaped domain, to be used in adaptive algorithms.
  */
  template <int d, int dt, int s0, int s1>
  class RingLaplacian
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
    
    //! constructor from a given wavelet basis and a given right-hand side y
    RingLaplacian(const WaveletBasis& basis,
		  const InfiniteVector<double,Index>& y);

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

  protected:
    const WaveletBasis& basis_;
    
    // rhs, mutable to have 'const' method
    mutable InfiniteVector<double,Index> y_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;

    // class for columns of 1D integral cache in angular direction
    typedef std::map<Index0,FixedArray1D<double,2> > Column1D_0;

    // class for columns of 1D integral cache in radial direction
    typedef std::map<Index1,FixedArray1D<double,4> > Column1D_1;

    // class for 1D integral cache in angular direction
    typedef std::map<Index0,Column1D_0> One_D_IntegralCache0;
    
    // class for 1D integral cache in radial direction
    typedef std::map<Index1,Column1D_1> One_D_IntegralCache1;

    //! cache for 1D integrals in angular direction
    mutable One_D_IntegralCache0 i12_cache;
    
    //! cache for 1D integrals in radial direction
    mutable One_D_IntegralCache1 i3456_cache;

    /*!
      helper function for RingLaplacian, compute various 1D integrals
      in angular direction (cf. RingLaplacian::a())
    */
    void angular_integrals(const Index0& lambda,
			   const Index0& mu,
			   double& I1, double& I2) const;

    /*!
      helper function for RingLaplacian, compute various 1D integrals
      in radial direction (cf. RingLaplacian::a())
    */
    void radial_integrals(const Index1& lambda,
			  const Index1& mu,
			  double& I3, double& I4, double& I5, double& I6) const;
  };
  
}

#include <galerkin/ring_laplacian.cpp>

#endif

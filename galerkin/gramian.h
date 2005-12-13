// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_GRAMIAN_H
#define _WAVELETTL_GRAMIAN_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <algebra/infinite_vector.h>

namespace WaveletTL
{
  /*!
    This class models equations Au=f with the biinfinite Gramian matrix 

      A = <Psi,Psi>

    of a given wavelet basis over a d-dimensional domain Omega,
    to be used in (adaptive) algorithms.

    The class has the minimal signature to be used within the APPLY routine
    or within adaptive solvers like CDD1.
  */
  template <class WBASIS, unsigned int DIM>
  class Gramian
  {
  public:
//     /*!
//       constructor from a given right-hand side
//     */
//     CubeEquation(const EllipticBVP<DIM>* bvp,
// 		 const FixedArray1D<bool,2*DIM>& bc);

//     /*!
//       copy constructor
//     */
//     CubeEquation(const CubeEquation&);

//     /*!
//       make template argument accessible
//     */
//     typedef CUBEBASIS WaveletBasis;

//     /*!
//       wavelet index class
//     */
//     typedef typename WaveletBasis::Index Index;

//     /*!
//       read access to the basis
//     */
//     const CUBEBASIS& basis() const { return basis_; }
    
//     /*!
//       space dimension of the problem
//     */
//     static const int space_dimension = DIM;

//     /*!
//       differential operators are local
//     */
//     static bool local_operator() { return true; }

//     /*!
//       (half) order t of the operator
//     */
//     static int operator_order() { return 1; }
    
//     /*!
//       evaluate the diagonal preconditioner D
//     */
//     double D(const typename WaveletBasis::Index& lambda) const;

//     /*!
//       rescale a coefficient vector by an integer power of D, c |-> D^{n}c
//     */
//     void rescale(InfiniteVector<double,typename WaveletBasis::Index>& coeffs,
// 		 const int n) const;

//     /*!
//       evaluate the (unpreconditioned) bilinear form a;
//       you can specify the order p of the quadrature rule, i.e.,
//       (piecewise) polynomials of maximal degree p will be integrated exactly.
//       Internally, we use an m-point composite tensor product Gauss rule adapted
//       to the singular supports of the spline wavelets involved,
//       so that m = (p+1)/2;
//     */
//     double a(const typename WaveletBasis::Index& lambda,
// 	     const typename WaveletBasis::Index& nu,
// 	     const unsigned int p = 4) const;

//     /*!
//       estimate the spectral norm ||A||
//     */
//     double norm_A() const;

//     /*!
//       estimate the spectral norm ||A^{-1}||
//     */
//     double norm_Ainv() const;

//     /*!
//       estimate compressibility exponent s^*
//       (we assume that the coefficients a(x),q(x) are smooth)
//     */
//     double s_star() const;
    
//     /*!
//       estimate the compression constants alpha_k in
//         ||A-A_k|| <= alpha_k * 2^{-s*k}
//     */
//     double alphak(const unsigned int k) const {
//       return 2*norm_A(); // suboptimal
//     }

//     /*!
//       evaluate the (unpreconditioned) right-hand side f
//     */
//     double f(const typename WaveletBasis::Index& lambda) const;

//     /*!
//       approximate the wavelet coefficient set of the preconditioned right-hand side F
//       within a prescribed \ell_2 error tolerance
//     */
//     void RHS(const double eta,
// 	     InfiniteVector<double,typename WaveletBasis::Index>& coeffs) const;

//     /*!
//       compute (or estimate) ||F||_2
//     */
//     double F_norm() const { return sqrt(fnorm_sqr); }

//     /*!
//       set the boundary value problem
//     */
//     void set_bvp(const EllipticBVP<DIM>*);

  protected:
//     const EllipticBVP<DIM>* bvp_;
//     CUBEBASIS basis_;

//     // right-hand side coefficients on a fine level, sorted by modulus
//     Array1D<std::pair<typename WaveletBasis::Index,double> > fcoeffs;

//     // precompute the right-hand side
//     void compute_rhs();

//     // (squared) \ell_2 norm of the precomputed right-hand side
//     double fnorm_sqr;

//     // estimates for ||A|| and ||A^{-1}||
//     mutable double normA, normAinv;
  };
}

#include <galerkin/gramian.cpp>

#endif

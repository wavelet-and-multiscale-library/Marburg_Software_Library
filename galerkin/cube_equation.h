// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_EQUATION_H
#define _WAVELETTL_CUBE_EQUATION_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>

#include <galerkin/galerkin_utils.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;

namespace WaveletTL
{
  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
      Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating an elliptic boundary value problem on the cube [0,1]^d
    
      -div(a(x)grad u(x)) + q(x)u(x) = f(x)

    with first (Dirichlet) or second (Neumann) order b.c.'s (modeled in
    the class EllipticBVP) as an equivalent operator equation
    within \ell_2 by means of a wavelet basis.

    The corresponding bilinear form in

      L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

      a(u,v) = \int_Omega [a(x)grad u(x)grad v(x)+q(x)u(x)v(x)] dx
      
    and the right-hand side is
 
      F(v) = \int_0^1 f(x)v(x) dx
   
    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from a wavelet basis \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.
    To achieve independence from the concrete choice of \Psi, the wavelet basis
    class is given as a template parameter CUBEBASIS. It should have a constructor of
    the form

       CUBEBASIS::CUBEBASIS(const FixedArray1D<bool,2*DIM>& bc);

    where bc indicates the enforcement of homogeneous Dirichlet boundary conditions (true).
    A natural concrete value for CUBEBASIS is the CubeBasis<DSBasis<d,dT> >.
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  class CubeEquation
  {
  public:
    /*!
      constructor from a boundary value problem and specified b.c.'s
    */
    CubeEquation(const EllipticBVP<DIM>* bvp,
		 const FixedArray1D<bool,2*DIM>& bc);

    /*!
      make template argument accessible
    */
    typedef CUBEBASIS WaveletBasis;

    /*!
      read access to the basis
    */
    const CUBEBASIS& basis() const { return basis_; }
    
    /*!
      space dimension of the problem
    */
    static int space_dimension() { return DIM; }

    /*!
      differential operators are local
    */
    static bool local_operator() { return true; }

    /*!
      (half) order t of the operator
    */
    static int operator_order() { return 1; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const typename WaveletBasis::Index& lambda) const;

    /*!
      rescale a coefficient vector by an integer power of D, c |-> D^{n}c
    */
    void rescale(InfiniteVector<double,typename WaveletBasis::Index>& coeffs,
		 const int n) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite tensor product Gauss rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
    double a(const typename WaveletBasis::Index& lambda,
	     const typename WaveletBasis::Index& nu,
	     const unsigned int p = 4) const;

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
      (we assume that the coefficients a(x),q(x) are smooth)
    */
    double s_star() const;
    
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
    double f(const typename WaveletBasis::Index& lambda) const;

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
	     InfiniteVector<double,typename WaveletBasis::Index>& coeffs) const;

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }

    /*!
      set the boundary value problem
    */
    void set_bvp(const EllipticBVP<DIM>*);

  protected:
    const EllipticBVP<DIM>* bvp_;
    CUBEBASIS basis_;

    // right-hand side coefficients on a fine level, sorted by modulus
    Array1D<std::pair<typename WaveletBasis::Index,double> > fcoeffs;

    // precompute the right-hand side
    void compute_rhs();

    // (squared) \ell_2 norm of the precomputed right-hand side
    double fnorm_sqr;

    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };
}

#include <galerkin/cube_equation.cpp>

#endif

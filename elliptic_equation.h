// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_ELLIPTIC_EQUATION_H
#define _FRAMETL_ELLIPTIC_EQUATION_H

#include <aggregated_frame.h>
#include <numerics/bvp.h>

using FrameTL::AggregatedFrame;
using MathTL::EllipticBVP;

namespace FrameTL
{

  /*!
    quadrature strategies for computation of
    stiffness matrix entries
  */
  enum QuadratureStrategy
    {
      Composite,
      TrivialAffine,
      SplineInterpolation
    };

  /*!
    This class models the (preconditioned) infinite-dimensional matrix problem
    
    Au = D^{-1}LD^{-1}u = D^{-1}F

    when reformulating a symmetric, second-order elliptic
    boundary value problem in divergence form over some domain
    Omega in R^d with boundary Gamma=dOmega,
    with homogeneous Dirichlet/Neumann/Robin boundary conditions

    -div(a(x)grad u(x)) + q(x)u(x) = f(x) in Omega
                             u(x) = 0 on Gamma_D
                         du/dn(x) = 0 on Gamma\Gamma_D.
    
    The corresponding bilinear form in

    L = (a(\psi_\nu,\psi_\lambda))_{\lambda,\nu}

    is

    a(u,v) = \int_Omega <a(x)*grad u(x), grad v(x)>  dx +
              \int_Omega q(x) * u(x) * v(x) dx
     
    and the right-hand side is
     
    f(v) = \int_Omega f(x)*v(x)  dx.

    The evaluation of a(.,.) and f is possible for arguments \psi_\lambda
    which stem from an aggregated wavelet frame \Psi=\{\psi_\lambda\} of the corresponding
    function space over Omega.     
  */
  template <class IBASIS, unsigned int DIM>
  class EllipticEquation
  {
  public:

//     /*!
//       constructor from a boundary value problem and specified b.c.'s
//     */
//     EllipticEquation(const EllipticBVP<DIM>* bvp,
// 		     const FixedArray1D<bool,2*DIM>& bc);

    /*!
      constructor
     */
    EllipticEquation(const EllipticBVP<DIM>* ell_bvp,
		     const AggregatedFrame<IBASIS,DIM>* frame,
		     const QuadratureStrategy qsrtat = Composite);

    /*!
      make template argument accessible
    */
    typedef AggregatedFrame<IBASIS,DIM> Frame;

    /*!
      make template argument accessible
    */
    typedef typename Frame::Index Index;

    /*!
      read access to the frame
    */
    const AggregatedFrame<IBASIS,DIM>*& frame() const { return fame; }
    
    /*!
      space dimension of the problem
    */
    static int space_dimension() { return DIM; }

    /*!
      differential operators are local
    */
    static bool local_operator() { return true; }

    /*!
      order of the operator
    */
    static int operator_order() { return 2; }
    
    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const;

    /*!
      rescale a coefficient vector by an integer power of D, c |-> D^{n}c
    */
    void rescale(InfiniteVector<double, typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs,
		 const int n) const;

    /*!
      evaluate the (unpreconditioned) bilinear form a;
      you can specify the order p of the quadrature rule, i.e.,
      (piecewise) polynomials of maximal degree p will be integrated exactly.
      Internally, we use an m-point composite Gauss quadrature rule adapted
      to the singular supports of the spline wavelets involved,
      so that m = (p+1)/2;
    */
    double a(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
	     const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
	     const unsigned int p = 2) const;

    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const;
    
    /*!
      returns spectral norm ||A^{-1}||
      estimate for ||A^{-1}|| has to be
      externally computed and to be set
      during initialization of the program.
      We assume this because ||A^{-1}|| quantity is hardly
      available in the frame case and
      quite complicated eigenvalue/eigenvector
      methods have to applied that are not implemented so
      far.
    */
    double norm_Ainv() const { return normAinv; };

    /*!
      sets estimate for ||A^{-1}||
    */
    double set_Ainv(const int nAinv) const { normAinv = nAinv; };

    /*!
      estimate compressibility exponent s^*
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
    double f(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda) const;

    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta, InfiniteVector<double, 
	     typename AggregatedFrame<IBASIS,DIM>::Index>& coeffs) const;

    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return sqrt(fnorm_sqr); }

    /*!
      set the boundary value problem
    */
    void set_bvp(const EllipticBVP<DIM>*);

   protected:
    
    /*!
      corresponding elliptic boundary value problem
     */
    const EllipticBVP<DIM>* ell_bvp_;

    /*!
      underlying frame
     */
    const AggregatedFrame<IBASIS,DIM>* frame_;

  private:

    /*!
      helper routines for a (.. , ..). Entries in diagonal and non-diagonal
      blocks of the stiffness matrix have to be treated differently.
     */
    double a_same_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			  const unsigned int q_order = 4) const;

    /*!
      helper routines for a (.. , ..). Entries in diagonal and non-diagonal
      blocks of the stiffness matrix have to be treated differently.
     */
    double a_different_patches(const typename AggregatedFrame<IBASIS,DIM>::Index& lambda,
			       const typename AggregatedFrame<IBASIS,DIM>::Index& nu,
			       const unsigned int q_order = 4) const;

    // precompute the right-hand side
    void compute_rhs();

    // right-hand side coefficients on a fine level, sorted by modulus
    Array1D<std::pair<typename AggregatedFrame<IBASIS,DIM>::Index,double> > fcoeffs;

    // (squared) \ell_2 norm of the precomputed right-hand side
    double fnorm_sqr;

    // reminder: This keyword can only be applied to non-static
    // and non-const data members of a class. If a data member is declared mutable,
    // then it is legal to assign a value to this data member from
    // a const member function.
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;

    QuadratureStrategy qstrat_; 

  };
}

#include <elliptic_equation.cpp>

#endif

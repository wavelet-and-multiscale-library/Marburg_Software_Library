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

    /*!
      constructor
     */
    EllipticEquation(const EllipticBVP<DIM>& ell_bvp,
		     const AggregatedFrame<IBASIS,DIM>* frame);
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
	     const unsigned int p = 4) const;

    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const;
    
    /*!
      estimate the spectral norm ||A^{-1}||
    */
    double norm_Ainv() const;

//     /*!
//       estimate compressibility exponent s^*
//     */
//     double s_star() const {
//       return 1.0 + WBASIS::primal_vanishing_moments(); // [St04a], Th. 2.3 for n=1
//     }

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
//     double f(const typename WBASIS::Index& lambda) const;

//     /*!
//       approximate the wavelet coefficient set of the preconditioned right-hand side F
//       within a prescribed \ell_2 error tolerance
//     */
//     void RHS(const double eta, InfiniteVector<double, typename WBASIS::Index>& coeffs) const;

//     /*!
//       compute (or estimate) ||F||_2
//     */
//     double F_norm() const { return sqrt(fnorm_sqr); }

//   protected:
//     const SimpleSturmBVP& bvp_;
//     WBASIS basis_;

//     // right-hand side coefficients on a fine level, sorted by modulus
//     Array1D<std::pair<typename WBASIS::Index,double> > fcoeffs;

//     // \ell_2 norm of the precomputed right-hand side
//     double fnorm_sqr;


  };
}

#include <elliptic_equation.cpp>

#endif

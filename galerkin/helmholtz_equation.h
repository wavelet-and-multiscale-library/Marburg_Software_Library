// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_HELMHOLTZ_EQUATION_H
#define _WAVELETTL_HELMHOLTZ_EQUATION_H

#include <utils/function.h>
#include <interval/spline_basis.h>

using MathTL::Function;

namespace WaveletTL
{
  /*!
    This class models the wavelet discretization of the one--dimensional
    Helmholtz equation

      -u''(t) + alpha * u(t) = f(t), 0 <= t <= 1

    with Dirichlet boundary conditions u(0)=u(1)=0.

    The wavelet bases useable for the discretization are those
    modeled by SplineBasis<d,dT>.
  */
  template <int d, int dT>
  class HelmholtzEquation1D
  {
  public:
    /*!
      constructor from a right-hand side f
    */
    HelmholtzEquation1D(const Function<1>* f,
			const double alpha = 0,
			const bool precompute_rhs = true);

    /*!
      type of the wavelet basis
    */
    typedef SplineBasis<d,dT> WaveletBasis;

    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return basis_; }
    
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
    
  protected:
    const Function<1>* f_;
    double alpha_;
    SplineBasis<d,dT> basis_;
  };


//     /*!
//       evaluate the diagonal preconditioner D
//     */
//     double D(const Index& lambda) const;

//     /*!
//       evaluate the (unpreconditioned) bilinear form a
//       (inherited from FullyDiagonalEnergyNormPreconditioner)
//     */
//     double a(const Index& lambda,
// 	     const Index& nu) const;

//     /*!
//       evaluate the (unpreconditioned) bilinear form a;
//       you can specify the order p of the quadrature rule, i.e.,
//       (piecewise) polynomials of maximal degree p will be integrated exactly.
//       Internally, we use an m-point composite Gauss quadrature rule adapted
//       to the singular supports of the spline wavelets involved,
//       so that m = (p+1)/2;
//     */
//     double a(const Index& lambda,
// 	     const Index& nu,
// 	     const unsigned int p) const;

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
//     double f(const Index& lambda) const;

//     /*!
//       approximate the wavelet coefficient set of the preconditioned right-hand side F
//       within a prescribed \ell_2 error tolerance
//     */
//     void RHS(const double eta, InfiniteVector<double,Index>& coeffs) const;

//     /*!
//       compute (or estimate) ||F||_2
//     */
//     double F_norm() const { return sqrt(fnorm_sqr); }

//   protected:
//     const SimpleSturmBVP& bvp_;
//     WBASIS basis_;

//     // flag whether right-hand side has already been precomputed
//     mutable bool rhs_precomputed;
    
//     /*!
//       precomputation of the right-hand side
//       (constness is not nice but necessary to have RHS a const function)
//     */
//     void precompute_rhs() const;

//     // right-hand side coefficients on a fine level, sorted by modulus
//     mutable Array1D<std::pair<Index,double> > fcoeffs;

//     // (squared) \ell_2 norm of the precomputed right-hand side
//     mutable double fnorm_sqr;

//     // estimates for ||A|| and ||A^{-1}||
//     mutable double normA, normAinv;
//   };
}

#include <galerkin/helmholtz_equation.cpp>

#endif

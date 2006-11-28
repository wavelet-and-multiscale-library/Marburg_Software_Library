// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_ROW_STAGE_EQUATION_H
#define _WAVELETTL_ROW_STAGE_EQUATION_H

#include <cmath>
#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <numerics/ivp.h>
#include <numerics/row_method.h>
#include <utils/function.h>
#include <galerkin/cached_problem.h>
#include <galerkin/infinite_preconditioner.h>
#include <adaptive/compression.h>
#include <galerkin/full_helmholtz.h>

using MathTL::Vector;
using MathTL::InfiniteVector;
using MathTL::AbstractIVP;
using MathTL::Function;
using MathTL::WMethodPreprocessRHSHelper;

namespace WaveletTL
{
  /*!
    This class models the operator and the right--hand side
    of the elliptic ROW stage equation

      D_alpha^{-1}<(alpha*I-A)Psi,Psi>D_alpha^{-1} * D_alpha*u_i
      = D_alpha^{-1}<r_{i,h},Psi>
      
    when discretizing a linear parabolic equation of the form

      u'(t) = Au(t) + f(t) =: F(t,u(t)),  0 < t <= T
       u(0) = u_0
       
    with a biorthogonal wavelet basis. Here A:H->H' is an isomorphism,
    and f: (0,T] -> H' is some driving term.
  
    The diagonally preconditioned stiffness matrix corresponding to A
      -A = D^{-1}(a(\psi_\lambda',\psi_\lambda))_{\lambda,\lambda'}D^{-1}
    is modeled in the template parameter class ELLIPTIC_EQ.
    
  */
  template <class ELLIPTIC_EQ, class GRAMIAN>
  class ROWStageEquation
    : public FullyDiagonalEnergyNormPreconditioner<typename ELLIPTIC_EQ::Index>
  {
  public:
    /*!
      Constructor from a given ROW method (finite-vector-valued, we only need the coeffs),
      a helper object for the stiffness matrix and the gramian,
      a (time-dependent) driving term f and its derivative f'.
      If ft and/or f are omitted (or set to zero), we assume that f and/or ft=0.
    */
    ROWStageEquation(const ROWMethod<Vector<double> >* row_method,
		     const ELLIPTIC_EQ* elliptic,
		     const GRAMIAN* gramian,
		     MathTL::Function<ELLIPTIC_EQ::space_dimension,double>* f = 0,
		     MathTL::Function<ELLIPTIC_EQ::space_dimension,double>* ft = 0)
      : row_method_(row_method),
	elliptic_(elliptic),
	gramian_(gramian),
	f_(f),
	ft_(ft),
	alpha_(0.0) // dummy value, alpha has to be set just before solving the stage eq.
    {
    }

    /*!
      make wavelet basis type accessible
    */
    typedef typename ELLIPTIC_EQ::WaveletBasis WaveletBasis;
    
    /*!
      wavelet index class
    */
    typedef typename ELLIPTIC_EQ::Index Index;
    
    /*!
      read access to the basis
    */
    const WaveletBasis& basis() const { return elliptic_->basis(); }
    
    /*!
      space dimension of the problem
    */
    static const int space_dimension = ELLIPTIC_EQ::space_dimension;
    
    /*!
      locality of the operator
    */
    static bool local_operator() { return ELLIPTIC_EQ::local_operator(); }
    
    /*!
      (half) order t of the operator
    */
    double operator_order() const { return elliptic_->operator_order(); }

    /*!
      evaluate the diagonal preconditioner D_alpha
    */
    double D(const Index& lambda) const { return sqrt(a(lambda, lambda)); }
    
    /*!
      evaluate the (unpreconditioned) bilinear form a
    */
    double a(const Index& lambda,
	     const Index& nu) const
    {
      return alpha_ * gramian_->a(lambda, nu) + elliptic_->a(lambda, nu);
    }
    
    /*!
      estimate the spectral norm ||D_alpha^{-1}*(alpha*G-A)*D_alpha^{-1}|| from above
    */
    double norm_A() const
    {
      return elliptic_->norm_A(); // reasonable for a large range of alpha
    }
      
    /*!
      estimate the spectral norm ||(D_alpha^{-1}*(alpha*G-A)*D_alpha^{-1})^{-1}|| from above
    */
    double norm_Ainv() const
    {
      return elliptic_->norm_Ainv(); // reasonable for a large range of alpha
    }

    /*!
      estimate compressibility exponent s^*
    */
    double s_star() const {
      return elliptic_->s_star(); // reasonable
    }
    
    /*!
      estimate the compression constants alpha_k in
      ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return elliptic_->alphak(k); // reasonable
    }
    
    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const {
      return y.get_coefficient(lambda);
    }
    
    /*!
      approximate the wavelet coefficient set of the preconditioned
      right-hand side F of the stage equation
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
	     InfiniteVector<double, Index>& coeffs) const {
      coeffs = y; // dirty
      coeffs.scale(this, -1); // coeffs *= D_alpha^{-1}
    }
    
    /*!
      compute (or estimate) ||F||_2
    */
    double F_norm() const { return l2_norm(y); }
    
    /*!
      w += factor * (stiffness matrix entries in column lambda on level j)
    */
    void add_level (const Index& lambda,
 		    InfiniteVector<double, Index>& w, const int j,
 		    const double factor,
 		    const int J,
 		    const CompressionStrategy strategy = St04a) const;
    
    /*!
      set (numerical) diffusion coefficient alpha
    */
    void set_alpha(double alpha) { alpha_ = alpha; }

    /*!
      setup (unpreconditioned) right-hand side y for the i-th stage equation,
      this requires that the stage solutions D_alpha u_1, ..., D_alpha u_{i-1}
      are already set up
    */
    void setup_rhs(unsigned int i,
		   const double tolerance,
		   const double t_n,
		   const double h,
		   const InfiniteVector<double,Index>& D_un,
		   const std::list<InfiniteVector<double,Index> >& Dalpha_uj,
		   const int jmax = 10);

    const ROWMethod<Vector<double> >* row_method_;
    const ELLIPTIC_EQ* elliptic_;
    const GRAMIAN* gramian_;
    
  protected:
    Function<space_dimension>* f_;
    Function<space_dimension>* ft_;

    double alpha_;

    //! holds current (unpreconditioned) right--hand side 
    InfiniteVector<double, typename ELLIPTIC_EQ::Index> y;
  };


  /*!
    an analogous class for approximations in full spaces V_j,
    using spline wavelets on the interval
  */
  template<int d, int dT>
  class ROWStageEquationFull1D
  {
  public:
    /*!
      Constructor from a given ROW method (finite-vector-valued, we only need the coeffs)
      and various helper objects,
      a (time-dependent) driving term f and its derivative f'.
      If ft and/or f are omitted (or set to zero), we assume that f and/or ft=0.
    */
    ROWStageEquationFull1D(const ROWMethod<Vector<double> >* row_method,
			   FullHelmholtz<d,dT>* A,
			   MathTL::Function<1,double>* f = 0,
			   MathTL::Function<1,double>* ft = 0)
      : row_method_(row_method),
	A_(A),
	f_(f),
	ft_(ft)
    {
    }
    
  protected:
    const ROWMethod<Vector<double> >* row_method_;
    FullHelmholtz<d,dT>* A_;
    Function<1>* f_;
    Function<1>* ft_;
  };

}

#include <parabolic/row_stage_equation.cpp>

#endif

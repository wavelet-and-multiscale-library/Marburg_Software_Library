// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LIN_PAR_EQUATION_H
#define _WAVELETTL_LIN_PAR_EQUATION_H

#include <algebra/infinite_vector.h>
#include <numerics/ivp.h>

using MathTL::InfiniteVector;
using MathTL::AbstractIVP;

namespace WaveletTL
{
  /*!
    ROW stage equation helper class for linear parabolic equations (see below).
    The helper class models the preconditioned stage equation

      D^{-1}(alpha*I-T)D^{-1}*Dx = D^{-1}y    

    The class has the minimal signature to be used within the APPLY routine
    or within adaptive solvers like CDD1.
  */
  template <class ELLIPTIC_EQ>
  class LinParEqROWStageEquationHelper
  {
  public:
    /*!
      constructor from alpha, T and y
    */
    LinParEqROWStageEquationHelper
    (const double alpha,
     const ELLIPTIC_EQ* T,
     const InfiniteVector<double, typename ELLIPTIC_EQ::Index>& y);
    
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
    const WaveletBasis& basis() const { return T->basis(); }
    
    /*!
      space dimension of the problem
    */
    static int space_dimension() { return ELLIPTIC_EQ::space_dimension(); }
    
    /*!
      locality of the operator
    */
    static bool local_operator() { return ELLIPTIC_EQ::local_operator(); }
    
    /*!
      (half) order t of the operator
    */
    static int operator_order() { return ELLIPTIC_EQ::operator_order(); }

    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const Index& lambda) const { return T->D(lambda); }
    
    /*!
      rescale a coefficient vector by an integer power of D, c |-> D^{n}c
    */
    void rescale(InfiniteVector<double,Index>& coeffs,
		 const int n) const {
      T->rescale(coeffs, n);
    }
    
    /*!
      evaluate the (unpreconditioned) bilinear form a
    */
    double a(const Index& lambda,
	     const Index& nu) const
    {
      return (lambda == nu ? alpha : 0.) + T->a(lambda, nu);
    }
    
    /*!
      estimate the spectral norm ||A||
    */
    double norm_A() const
    {
      return T->norm_A(); // dirty
    }
      
    /*!
      estimate the spectral norm ||A^{-1}||
    */
    double norm_Ainv() const
    {
      return T->norm_Ainv(); // dirty
    }
    
    /*!
      estimate compressibility exponent s^*
    */
    double s_star() const {
      return T->s_star();
    }
    
    /*!
      estimate the compression constants alpha_k in
      ||A-A_k|| <= alpha_k * 2^{-s*k}
    */
    double alphak(const unsigned int k) const {
      return 2*norm_A(); // pessimistic
    }
    
    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const {
      return y.get_coefficient(lambda);
    }
    
    /*!
      approximate the wavelet coefficient set of the preconditioned right-hand side F
      within a prescribed \ell_2 error tolerance
    */
    void RHS(const double eta,
	     InfiniteVector<double, Index>& coeffs) const {
      coeffs = y; // dirty
      T->rescale(coeffs, -1);
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

  protected:
    const double alpha;
    const ELLIPTIC_EQ* T;
    const InfiniteVector<double, typename ELLIPTIC_EQ::Index> y;
  };

  /*!
    This class models a linear parabolic equation of the form

      u'(t) = Au(t) + f(t) =: F(t,u(t)),  0 < t <= T
       u(0) = u_0

    where -A: \ell_{2,D} -> \ell_{2,D^{-1}} is a positive isomorphism
    and f: (0,T] -> \ell_{2,D^{-1}}, D being a diagonal preconditioner (see below).
    An equation of this form arises when equivalently reformulating
    a problem of the analogous form

      v'(t) = Bv(t) + g(t) =: G(t,v(t)), 0 < t <= T
       v(0) = v_0

    with isomorphism -B: H -> H' and right-hand side g: (0,T] -> H'
    by means of a biorthogonal wavelet basis Psi=\{\psi_\lambda\},
    setting F(t,v):=<G(t,v^\top\Psi),\tilde\Psi> for all v in \ell_{2,D}.
    
    The (unpreconditioned) stiffness matrix

      -A = (a(\psi_\lambda',\psi_\lambda))_{\lambda,\lambda'}

    is modeled in the template parameter class ELLIPTIC_EQ.
    
  */
  template <class ELLIPTIC_EQ>
  class LinearParabolicEquation
    : public AbstractIVP<InfiniteVector<double, typename ELLIPTIC_EQ::Index> >
  {
  public:
    /*!
      the index class, for convenience
    */
    typedef typename ELLIPTIC_EQ::Index Index;

    /*!
      constructor from a helper object for the stiffness matrix
      and a given initial value u0 in ell_2
     */
    LinearParabolicEquation(const ELLIPTIC_EQ* elliptic,
			    const InfiniteVector<double,Index>& u0);

    /*!
      evaluate the right-hand side F(t,v)=Av+f(t) up to a prescribed tolerance
    */
    void evaluate_f(const double t,
		    const InfiniteVector<double,Index>& v,
		    const double tolerance,
		    InfiniteVector<double,Index>& result) const;
    
    /*!
      evaluate F_t(t,v)=f'(t) up to a prescribed tolerance
    */
    void evaluate_ft(const double t,
		     const InfiniteVector<double,Index>& v,
		     const double tolerance,
		     InfiniteVector<double,Index>& result) const;
    
    /*!
      solve (alpha*I-A)x=y up to a prescribed tolerance
    */
    void solve_ROW_stage_equation(const double t,
				  const InfiniteVector<double,Index>& v,
				  const double alpha,
				  const InfiniteVector<double,Index>& y,
				  const double tolerance,
				  InfiniteVector<double,Index>& result) const;

  protected:
    //! pointer to the elliptic subproblem helper
    const ELLIPTIC_EQ* elliptic;
  };
}

#include <parabolic/lin_par_equation.cpp>

#endif

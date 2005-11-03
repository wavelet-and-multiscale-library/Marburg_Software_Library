// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_W_METHOD_H
#define _MATHTL_W_METHOD_H

#include <algebra/triangular_matrix.h>
#include <algebra/vector.h>

namespace MathTL
{
  /*!
    The following class is an abstract base for an s-stage linearly implicit
    method (a W-method) for the numerical solution of
    (abstract) nonautonomous initial value problems of the form

      u'(t) = F(t, u(t)), u(0) = u_0.

    The numerical method solves the following s linear stage equations:

      (I-\tau*\gamma_{i,i}*T)(k_i + \sum_{j=1}^{i-1}\frac{\gamma_{i,j}}{\gamma_{i,i}}k_j)
        = F(t_m+\tau*\alpha_i, u^{(m)}+\tau * \sum_{j=1}^{i-1}\alpha_{i,j} * k_j)
	  + \sum_{j=1}^{i-1} \frac{\gamma_{i,j}}{\gamma_{i,i}} * k_j
	  + \tau * \gamma_i * g
	  
    Here \alpha_i = \sum_{j=1}^i\alpha_{i,j}, \gamma_i = \sum_{j=1}^i\gamma_{i,j}.
    T is the exact Jacobian F_u(t_m,u^{(m)}) or an approximation of it.
    In the latter case, one speaks of W-methods.
    g is the exact derivative F_t(t_m,u^{(m)}) or an approximation of it.
    In a general W-method, one will often choose g=0.
    The value at the new time node t_{m+1}=t_m+\tau is then determined by

      u^{(m+1)} = u^{(m)} + \tau * \sum_{i=1}^s b_i*k_i

    The template parameter VECTOR stands for an element of the Hilbert/Banach space X
    under consideration. For finite-dimensional problems, this will almost always
    be a vector class like Vector<double>. In the infinite-dimensional case,
    VECTOR models elements of X as (finite) linear combinations of a wavelet basis.

    All methods under consideration provide embedded lower order error estimators.
  */
  template <class VECTOR>
  class WMethod
    : public OneStepScheme<VECTOR>
  {
  public:
    /*!
      constructor from the W-method coefficients
    */
    WMethod();

    /*!
      virtual destructor
    */
    virtual ~WMethod() {}

    /*!
      increment function u^{(m)} -> u^{(m+1)},
      also returns a local error estimator
    */
    void increment(const AbstractIVP<VECTOR>& ivp,
		   const double t_m,
		   const VECTOR& u_m,
		   const double tau,
		   VECTOR& u_mplus1,
		   VECTOR& error_estimate,
		   const double tolerance = 1e-2) const;

    /*!
      (adaptive) solver for one of the systems (I-\alpha*T)x=y,
      implicitly this also specifies the approximation T to the Jacobian F_v(t_m,u^{(m)})
    */
    virtual void solve_stage_equation(const AbstractIVP<VECTOR>& ivp,
				      const double alpha,
				      const VECTOR& y,
				      const double tolerance,
				      VECTOR& x) const = 0;
    
    /*!
      evaluate the vector g (an approximation of F_t(t_m,u^{(m)}))
    */
    virtual void g(const AbstractIVP<VECTOR>& ivp,
		   const double t_m,
		   const VECTOR& u_m,
		   const double tolerance,
		   VECTOR& result) const = 0;
    
  protected:
    //! W-method coefficients
    LowerTriangularMatrix<double> alpha, gamma_matrix;
    Vector<double> b, gamma_vector;
  };
}

#include <numerics/w_method.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ROW_METHOD_H
#define _MATHTL_ROW_METHOD_H

#include <algebra/triangular_matrix.h>
#include <algebra/vector.h>
#include <numerics/w_method.h>

namespace MathTL
{
  /*!
    The following class is an abstract base for an s-stage Rosenbrock-Wanner method
    (ROW-method) for the numerical solution of (abstract)
    nonautonomous initial value problems of the form

      u'(t) = F(t, u(t)), u(0) = u_0.

    Essentially, a ROW-method is a W-method with T = F_v(t_m,u^{(m)}) and
    g = F_t(t_m,u^{(m)}. Please read w_method.h for details.
  */
  template <class VECTOR>
  class ROWMethod
    : public WMethod<VECTOR>
  {
  public:
    /*!
      constructor for one of the builtin ROW-methods, cf. w_method.h for details
    */
    ROWMethod(const typename WMethod<VECTOR>::Method method,
	      const WMethodStageEquationSolver<VECTOR>& stage_equation_solver);

    /*!
      virtual destructor
    */
    virtual ~ROWMethod() {}

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
      (adaptive) solver for one of the systems (I-\alpha*T)x=y
    */
    void solve_stage_equation(const AbstractIVP<VECTOR>& ivp,
			      const double alpha,
			      const VECTOR& y,
			      const double tolerance,
			      VECTOR& x) const;
    
    /*!
      evaluate the vector g (an approximation of F_t(t_m,u^{(m)}))
    */
    void g(const AbstractIVP<VECTOR>& ivp,
	   const double t_m,
	   const VECTOR& u_m,
	   const double tolerance,
	   VECTOR& result) const;
  };
}

#include <numerics/row_method.cpp>

#endif

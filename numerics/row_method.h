// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ROW_METHOD_H
#define _MATHTL_ROW_METHOD_H

#include <algebra/triangular_matrix.h>
#include <algebra/vector.h>
#include <numerics/w_method.h>

namespace MathTL
{
  /*!
    The following class is an abstract base for an s-stage
    Rosenbrock-Wanner (ROW-) method
    for the numerical solution of (abstract)
    nonautonomous initial value problems of the form

      u'(t) = F(t, u(t)), u(0) = u_0.

    Essentially, a ROW-method is a W-method with T = F_v(t_m,u^{(m)}) and
    g = F_t(t_m,u^{(m)}, see w_method.h for details.
  */
  template <class VECTOR>
  class ROWMethod
    : public WMethod<VECTOR>, public WMethodStageEquationHelper<VECTOR>
  {
  public:
    /*!
      constructor for one of the builtin ROW-methods, cf. w_method.h for details
    */
    ROWMethod(const typename WMethod<VECTOR>::Method method);

    /*!
      virtual destructor
    */
    virtual ~ROWMethod() {}

    /*!
      (adaptive) solver for one of the systems (alpha*I-T)x=y,
      inherited from WMethodStageEquationHelper
    */
    void solve_W_stage_equation(const AbstractIVP<VECTOR>* ivp,
				const double t,
				const VECTOR& v,
				const double alpha,
				const VECTOR& y,
				const double tolerance,
				VECTOR& x) const
    {
      // use the exact jacobian
      ivp->solve_ROW_stage_equation(t, v, alpha, y, tolerance, x);
    }

    /*!
      Approximation of the temporal derivative f_t(t_m,u^{(m)}),
      inherited from WMethodStageEquationHelper
    */
    void approximate_ft(const AbstractIVP<VECTOR>* ivp,
			const double t,
			const VECTOR& v,
			const double tolerance,
			VECTOR& result) const
    {
      // use the exact derivative f_t
      ivp->evaluate_ft(t, v, tolerance, result);
    }
  };
}

#include <numerics/row_method.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_IVP_H
#define _MATHTL_IVP_H

#include <geometry/point.h>

namespace MathTL
{
  /*!
    abstract base class for a vector-valued initial value problem
    
      u'(t) = f(t, u(t)),   0 < t <= T
      u(0) = u_0

    where u:[0,T]->\mathbb R^d

    Note that the class fulfills the necessary signature to be used in
    one of the Rosenbrock methods for the numerical approximation of u(t).
  */
  template <unsigned int DIM>
  class IVP
  {
  public:
    /*!
      initial value
    */
    Point<DIM> u0;

    /*!
      virtual destructor
     */
    virtual ~IVP ();
    
    /*!
      evaluate the right--hand side f at (t,v)
    */
    virtual void apply_f(const double t, const Point<DIM>& v,
			 Point<DIM>& result) const = 0;

    /*!
      evaluate the partial derivative f_t at (t,v)
    */
    virtual void apply_ft(const double t, const Point<DIM>& v,
			  Point<DIM>& result) const = 0;

    /*!
      solve the special linear system
      
        (alpha*I-J)u = v,
      
      where J = \partial_v f(t,v) is the Jacobian of f
    */
    virtual void solve_jacobian(const double t, const Point<DIM>& v, const double alpha,
				Point<DIM>& result) const = 0;
  };

  /*!
    Abstract base class for general initial value problems
    
      u'(t) = f(t, u(t)),   0 < t <= T
      u(0) = u_0

    where u:[0,T]->V.
    
    The signature of AbstractIVP is designed to be used in (derivations of) the
    class OneStepScheme, especially for Runge-Kutta and linearly implicit methods.
  */
  template <class VECTOR>
  class AbstractIVP
  {
  public:
    /*!
      initial value
    */
    VECTOR u0;

    /*!
      virtual destructor
    */
    virtual ~AbstractIVP() = 0;

    /*!
      evaluate the right--hand side f at (t,v),
      up to some tolerance (w.r.t. the ||.||_2 norm)
    */
    virtual void evaluate_f(const double t,
			    const VECTOR& v,
			    const double tolerance,
			    VECTOR& result) const = 0;

    /*!
      evaluate the derivative f_t at (t,v),
      up to some tolerance (w.r.t. the ||.||_2 norm)
    */
    virtual void evaluate_ft(const double t,
			     const VECTOR& v,
			     const double tolerance,
			     VECTOR& result) const = 0;

    /*!
      Up to a given tolerance (w.r.t. the ||.||_2 norm), solve the special
      linear system 
      
        (alpha*I-J)x = y,
      
      where J = \partial_v f(t,v) is the (exact) Jacobian of f.
      This routine will primarily be used to solve the stage equations of
      an ROW-method.
    */
    virtual void solve_ROW_stage_equation(const double t,
					  const VECTOR& v,
					  const double alpha,
					  const VECTOR& y,
					  const double tolerancs,
					  VECTOR& result) const = 0;
  };
}

#include <numerics/ivp.cpp>

#endif

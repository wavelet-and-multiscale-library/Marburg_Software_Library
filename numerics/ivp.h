// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
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
      apply the right--hand side f to (t,v)
    */
    virtual void apply_f(const double t, const Point<DIM>& v,
			 Point<DIM>& result) const = 0;

    /*!
      apply the partial derivative f_t to (t,v)
    */
    virtual void apply_ft(const double t, const Point<DIM>& v,
			  Point<DIM>& result) const = 0;

    /*!
      solve the special linear system
      
        (I-\tau*J)u = v,
      
      where J = \partial_v f(t,v) is the Jacobian of f
    */
    virtual void solve_jacobian(const double t, const Point<DIM>& v, const double tau,
				Point<DIM>& result) const = 0;
  };

  /*!
    Abstract base class for general initial value problems
    
      u'(t) = f(t, u(t)),   0 < t <= T
      u(0) = u_0

    where u:[0,T]->V.
    
    The signature of AbstractIVP is designed to be used in (derivations of) the
    class OneStepScheme.
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
      apply the right--hand side f to (t,v),
      up to some tolerance (in the ||.||_2 norm)
    */
    virtual void apply_f(const double t, const VECTOR& v, const double tolerance,
			 VECTOR& result) const = 0;
  };
}

#include <numerics/ivp.cpp>

#endif

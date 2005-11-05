// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ONESTEPSCHEME_H
#define _MATHTL_ONESTEPSCHEME_H

#include <numerics/ivp.h>

namespace MathTL
{
  /*!
    Abstract base class for various one-step schemes for the
    adaptive numerical solution of vector-valued initial value problems
    
      u'(t) = f(t, u(t)),   0 < t <= T
      u(0) = u_0

    where u:[0,T]->V. A one-step schemes essentially has to perform
    the increment u^{(m)}->u^{(m+1)}, where u^{(m)} is an approximation
    of u(t_m) for a given time sequence 0=t_0<t_1<...
    Additionally, the increment function has to provide some local error
    estimator.

    Since the one-step scheme may also be used in an infinite-dimensional
    setting, it is possible to specify a tolerance parameter (which controls
    the overall spatial error introduced by the evaluations of the
    right-hand side).
  */
  template <class VECTOR>
  class OneStepScheme
  {
  public:
    /*!
      virtual destructor
    */
    virtual ~OneStepScheme() = 0;

    /*!
      increment function u^{(m)} -> u^{(m+1)},
      also returns a local error estimator
    */
    virtual void increment(const AbstractIVP<VECTOR>* ivp,
			   const double t_m, const VECTOR& u_m,
			   const double tau,
			   VECTOR& u_mplus1,
			   VECTOR& error_estimate,
			   const double tolerance = 1e-2) const = 0;
  };
}

#include <numerics/one_step_scheme.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ONESTEPSCHEME_H
#define _MATHTL_ONESTEPSCHEME_H

#include <list>
#include <map>
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

    /*!
      consistency/convergence order of the scheme
    */
    virtual int order() const = 0;
  };

  /*!
    type of the approximate solution returned by the adaptive IVP solver
  */
  template <class VECTOR>
  class IVPSolution
  {
  public:
    std::list<double> t;
    std::list<VECTOR> u;
  };

  /*!
    Solve a given initial value problem on [0,T] adaptively with a given one-step scheme.
    You have to specify a maximal stepwidth tau_max and a factor q which limits the
    increase of the stepwidth when the error estimator is too small.
    The algorithm returns the time step sequence and the solution.

    Both absolute and relative tolerances can be specified, a step is accepted when
      ||y-yhat|| <= atol + max(||u_m||,||u_{m+1}||) * rtol
  */
  template <class VECTOR>
  void solve_IVP(const AbstractIVP<VECTOR>* ivp,
		 const OneStepScheme<VECTOR>* scheme,
		 const double T,
		 const double atol,
		 const double rtol,
		 const double q,
		 const double tau_max,
		 IVPSolution<VECTOR>& result);
}

#include <numerics/one_step_scheme.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_RUNGE_KUTTA_H
#define _MATHTL_RUNGE_KUTTA_H

#include <algebra/vector.h>
#include <algebra/triangular_matrix.h>
#include <numerics/one_step_scheme.h>

namespace MathTL
{
  /*!
    The following class models an s-stage explicit (embedded) Runge-Kutta method
    for the numerical solution of (abstract) nonautonomous initial value
    problems of the form

      u'(t) = F(t, u(t)), u(0) = u_0.
  */
  template <class VECTOR>
  class ExplicitRungeKuttaScheme
    : public OneStepScheme<VECTOR>
  {
  public:
    /*!
      enum type for different builtin RK schemes
      
      References:
      [DB] Deuflhard/Bornemann, Numerische Mathematik II, de Gruyter
      [S]  B.A.Schmitt, Numerik IIb, Vorlesungsskript
     */
    enum Method {
      RK12,       // s=2, p=2(1) (embedded scheme is Euler scheme)
      RK23,       // s=3, p=3(2) (embedded scheme is the Runge/improved Euler scheme), see [S]
      Fehlberg34, // s=5 (effective: 4), p=4(3), FSAL (uses the classical RK scheme), see [DB]
      DoPri45,    // s=7 (effective: 6), p=5(4), FSAL
      DoPri78     // s=13 (effective: 12), p=8(7), FSAL
    };

    /*!
      constructor from one of the builtin RK schemes
    */
    ExplicitRungeKuttaScheme(const Method method);

    /*!
      increment function + local error estimation
    */
    void increment(const AbstractIVP<VECTOR>* ivp,
		   const double t_m, const VECTOR& u_m,
		   const double tau,
		   VECTOR& u_mplus1,
		   VECTOR& error_estimate,
		   const double tolerance = 1e-2) const;

    /*!
      consistency/convergence order p
    */
    int order() const { return p; }

  protected:
    /*!
      "first same as last" flag
      (i.e., the first f-evaluation in the next time step is the last one of the current step)
    */
    bool fsal;

    //! node vector
    Vector<double> c;

    //! Butcher coefficients
    LowerTriangularMatrix<double> A;

    //! weight vectors
    Vector<double> b, bhat;

    //! consistency
    int p;
  };
}

#include "numerics/runge_kutta.cpp"

#endif

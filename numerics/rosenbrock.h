// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ROSENBROCK_H
#define _MATHTL_ROSENBROCK_H

#include <algebra/triangular_matrix.h>
#include <algebra/vector.h>

namespace MathTL
{
  /*!
    The following class models an s-stage Rosenbrock-type method
    for the numerical solution of (abstract) nonautonomous initial value
    problems of the form

      u'(t) = F(t, u(t)), u(0) = u_0.

    The numerical method solves the following s linear stage equations:

      (I-\tau*\gamma_{i,i}*J)(k_i + \sum_{j=1}^{i-1}\frac{\gamma_{i,j}}{\gamma_{i,i}}k_j)
        = F(t_m+\tau*\alpha_i, u^{(m)}+\tau * \sum_{j=1}^{i-1}\alpha_{i,j} * k_j)
	  + \sum_{j=1}^{i-1} \frac{\gamma_{i,j}}{\gamma_{i,i}} * k_j
	  + \tau * \gamma_i * \partial_tF(t_m,u^{(m)})
	  
    Here \alpha_i = \sum_{j=1}^i\alpha_{i,j}, \gamma_i = \sum_{j=1}^i\gamma_{i,j}.
    J is the exact Jacobian F_u(t_n, u_n) or an approximation of it.
    In the latter case, one speaks of W-methods.
    The value at the new time node t_{m+1}=t_m+\tau is then determined by

      u^{(m+1)} = u^{(m)} + \tau * \sum_{i=1}^s b_i*k_i

    The template parameter VECTOR stands for an element of the Banach space X
    under consideration. For finite-dimensional problems, this will almost always
    be a vector class like Vector<double>. In the infinite-dimensional case,
    VECTOR models elements of X as (finite) linear combinations of a wavelet basis.
    IVP provides the underlying initial value problem, this class should follow
    the signature

    class IVP
    {
      void apply_f(const double t, const VECTOR& v, VECTOR& result) const;
      void apply_ft(const double t, const VECTOR& v, VECTOR& result) const;
      void solve_jacobian(const double t, const VECTOR& v, const double tau,
			  VECTOR& result) const;
    }

    The methods under consideration provide embedded lower order error estimators.
  */
  template <class VECTOR, class IVP>
  class Rosenbrock
  {
  public:
    /*!
      enum type for the different builtin methods (in ascending stage order s)

      References:
      [Euler]  Deuflhard/Bornemann, Numerik II
      [ROS2]   Blom, Hundsdorfer, Spee, Verwer:
               A Second-Order Rosenbrock Method Applied to Photochemical Dispersion Problems,
	       SIAM J. Sci. Comput. 20(1999), 1456-1480
      [ROS3],
      [RODAS3] Blom, Carmichael, Potra, Sandu, Spee, Verwer:
               Benchmarking Stiff ODE Solvers for Atmospheric Chemistry Problems II:
	       Rosenbrock Solvers,
	       Atmos. Environ. 31(1997), 3459-3472
      [ROWDA3] Roche:
               Rosenbrock Methods for Differential Algebraic Equations,
               Numer. Math. 52(1988), 45-63
     */
    enum Method {
      Euler,  // s=1, p=1, L-stable
      ROS2,   // s=2, p=2, L-stable
      ROS3,   // s=3, p=3, L-stable, stiffly accurate
      ROWDA3, // s=3, p=3, L-stable
      RODAS3  // s=4, p=3, L-stable, stiffly accurate
    };

    /*!
      construct one of the builtin Rosenbrock methods
     */
    Rosenbrock(const Method method = Euler);

    /*!
      increment function u^{(m)} -> u^{(m+1)}
    */
    void increment(const IVP& ivp,
		   const double t_m, const VECTOR& u_m,
		   const double tau,
		   VECTOR& u_mplus1) const;
    
    /*!
      increment function u^{(m)} -> u^{(m+1)},
      also returns a local error estimator
    */
    void increment(const IVP& ivp,
		   const double t_m, const VECTOR& u_m,
		   const double tau,
		   VECTOR& u_mplus1,
		   VECTOR& error_estimate) const;
  protected:
    /*!
      the Rosenbrock coefficients
    */
    LowerTriangularMatrix<double> alpha_, gamma_;
    Vector<double> b_;

    /*!
      toggle method
    */
    Method method_;
  };
}

// include implementation
#include "numerics/rosenbrock.cpp"

#endif

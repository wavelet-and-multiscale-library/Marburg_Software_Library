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
    Abstract helper class for (adaptive) solvers of the W-(or ROW-) method
    stage equations

      (I-alpha*T)x=y.
      
    In a W-method, one replaces the exact Jacobian J=f_v(t,v) by an maybe crude
    approximation T. The derivative f_t(t,v) will also be approximated by a vector
    g (often even by zero).

    In practice, it is possible (but not necessary) to make the concrete instance
    of AbstractIVP under consideration also a descendant from this helper class.
    For example, the class ROWMethod uses a special WMethodStageEquationHelper
    class, see row_method.h for details.
  */
  template <class VECTOR>
  class WMethodStageEquationHelper
  {
  public:
    //! purely virtual destructor
    virtual ~WMethodStageEquationHelper() = 0;
    
    /*!
      (adaptive) solver for one of the systems (I-\alpha*T)x=y,
      implicitly this also specifies the approximation T to the
      Jacobian f_v(t,v)
    */
    virtual void solve_W_stage_equation(const AbstractIVP<VECTOR>* ivp,
					const double t,
					const VECTOR& v,
					const double alpha,
					const VECTOR& y,
					const double tolerance,
					VECTOR& x) const = 0;

    /*!
      Approximation of the temporal derivative f_t(t,v)
      (up to some tolerance w.r.t. ||.||_2, often denoted as g).
      In practice, g will often be chosen as zero, for stability reasons
      (reducing the overall consistency order of the W-method).
    */
    virtual void approximate_ft(const AbstractIVP<VECTOR>* ivp,
				const double t,
				const VECTOR& v,
				const double tolerance,
				VECTOR& result) const = 0;
  };
  
  /*!
    The following class is an abstract base for an s-stage linearly implicit
    method (a W-method) for the numerical solution of
    (abstract) nonautonomous initial value problems of the form

      u'(t) = f(t, u(t)), u(0) = u_0.

    The numerical method solves the following s linear stage equations:

      (I-\tau*\gamma_{i,i}*T)(k_i + \sum_{j=1}^{i-1}\frac{\gamma_{i,j}}{\gamma_{i,i}}k_j)
        = f(t_m + \tau * \alpha_i, u^{(m)} + \tau * \sum_{j=1}^{i-1}\alpha_{i,j} * k_j)
	  + \sum_{j=1}^{i-1} \frac{\gamma_{i,j}}{\gamma_{i,i}} * k_j
	  + \tau * \gamma_i * g
	  
    Here \alpha_i = \sum_{j=1}^{i-1}\alpha_{i,j}, \gamma_i = \sum_{j=1}^i\gamma_{i,j}.
    T is the exact Jacobian f_v(t_m,u^{(m)}) (ROW-method) or an approximation of it (W-method).
    g is the exact derivative f_t(t_m,u^{(m)}) (ROW-method) or an approximation of it (W-method).
    In a general W-method, one will often choose g=0.
    The value at the new time node t_{m+1}=t_m+\tau is then determined by

      u^{(m+1)} = u^{(m)} + \tau * \sum_{i=1}^s b_i*k_i

    The template parameter VECTOR stands for an element of the Hilbert/Banach space X
    under consideration. For finite-dimensional problems, this will almost always
    be a vector class like Vector<double>. In the infinite-dimensional case,
    VECTOR models elements of X as (finite) linear combinations of a wavelet basis.

    All methods under consideration provide embedded lower order error estimators,
    of the form

      uhat^{(m+1)} = u^{(m)} + \tau * \sum_{i=1}^s bhat_i*k_i
  */
  template <class VECTOR>
  class WMethod
    : public OneStepScheme<VECTOR>
  {
  public:
    /*!
      enum type for the different builtin methods

      References:
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
      ROS2,      // s=2, p=2, L-stable
      ROWDA3,    // s=3, p=3, L-stable
      ROS3,      // s=3, p=3, L-stable
      RODAS3     // s=4, p=3, L-stable, stiffly accurate
    };

    /*!
      constructor from some builtin methods and a given stage equation solver
    */
    WMethod(const Method method,
	    const WMethodStageEquationHelper<VECTOR>* stage_equation_solver);
    
    /*!
      virtual destructor
    */
    virtual ~WMethod() {}

    /*!
      increment function u^{(m)} -> u^{(m+1)},
      also returns a local error estimator
    */
    void increment(const AbstractIVP<VECTOR>* ivp,
		   const double t_m,
		   const VECTOR& u_m,
		   const double tau,
		   VECTOR& u_mplus1,
		   VECTOR& error_estimate,
		   const double tolerance = 1e-2) const;

    /*!
      consistency/convergence order
    */
    int order() const { return p; }
    
  protected:
    //! (alpha_{i,j})
    LowerTriangularMatrix<double> alpha_matrix;

    //! (gamma_{i,j})
    LowerTriangularMatrix<double> gamma_matrix;

    //! (b_i)
    Vector<double> b;

    //! (alpha_i), for convenience
    Vector<double> alpha_vector;

    //! (gamma_i), for convenience
    Vector<double> gamma_vector;

    //! coefficients for the lower order embedded scheme (bhat_i)
    Vector<double> bhat;

    //! instance of the stage equation helper
    const WMethodStageEquationHelper<VECTOR>* stage_equation_helper;

    //! consistency order
    int p;
  };
}

#include <numerics/w_method.cpp>

#endif

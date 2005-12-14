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
#include <numerics/ivp.h>
#include <numerics/one_step_scheme.h>

namespace MathTL
{
  /*!
    Abstract helper class for (adaptive) solvers of the W-(or ROW-) method
    stage equations

      (alpha*I-T)x=y.
      
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
      (adaptive) solver for one of the systems (alpha*I-T)x=y,
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
    Helper class for preprocessing the vector
      w:=\sum_{j=1}^{i-1} \frac{c_{i,j}}{\tau} * u_j
    by, e.g., a linear transformation.
    Nontrivial instances of this class are needed when it comes to
    a discretization of PDEs with biorthogonal wavelet bases. There the
    vector w has to be multiplied with the Gramian <Psi,Psi>.
  */
  template <class VECTOR>
  class WMethodPreprocessRHSHelper
  {
  public:
    //! purely virtual destructor
    virtual ~WMethodPreprocessRHSHelper() = 0;

    //! preprocess the vector w (up to the current tolerance)
    virtual void preprocess_rhs_share(VECTOR& wbeforeandafter,
				      const double tolerance) const
    { // default behaviour: do nothing
    }
  };
  
  /*!
    The following class is an abstract base for an s-stage linearly implicit
    method (a W-method) for the numerical solution of
    (abstract) nonautonomous initial value problems of the form

      u'(t) = f(t, u(t)), u(0) = u_0.

    The numerical method solves the following s linear stage equations:

      ((\tau*\gamma_{i,i})^{-1}I - T) u_i
        = f(t_m + \tau * \alpha_i, u^{(m)} + \sum_{j=1}^{i-1} a_{i,j} * u_j)
	  + \sum_{j=1}^{i-1} \frac{c_{i,j}}{\tau} * u_j
	  + \tau * \gamma_i * g
	  
    This form arises when reformulating the original stage equations

      (I-\tau*\gamma_{i,i}*T)k_i
        = f(t_m + \tau * \alpha_i, u^{(m)} + \tau * \sum_{j=1}^{i-1} \alpha_{i,j} * k_j)
	  + \tau * T\sum_{j=1}^{i-1} \gamma_{i,j} * k_j
	  + \tau * \gamma_i * g

    by introducing the new variables u_i=\tau\sum_{j=1}^{i-1}\gamma_{i,j}k_j,
    C=(c_{i,j})=diag(Gamma)^{-1}-Gamma^{-1}, A=(a_{i,j})=(alpha_{i,j})Gamma^{-1}.
    After this transformation, quite a few multiplications can be saved.
    The original stages and the new solution can be recovered via

      k_i = 1/gamma_{i,i}*u_i - \sum_{j=1}^{i-1}c_{i,j}u_j

      u^{(m+1)} = u^{(m)} + \tau * \sum_{i=1}^s b_i*k_i
                = u^{(m)} + \sum_{i=1}^s m_i*u_i,

    where (m_1,...,m_s)=(b_1,...,b_s)*Gamma^{-1}.
    The values gamma_{i,i} are stored in the diagonal of C.
    We assume that \alpha_i = \sum_{j=1}^{i-1}\alpha_{i,j}, \gamma_i = \sum_{j=1}^i\gamma_{i,j}.
    T is the exact Jacobian f_v(t_m,u^{(m)}) (ROW-method) or an approximation of it (W-method).
    g is the exact derivative f_t(t_m,u^{(m)}) (ROW-method) or an approximation of it (W-method).
    In a general W-method, one will often choose g=0.

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
      [GRK4T]  Kaps, Rentrop:
               Generalized Runge-Kutta methods of order four with stepsize control
	       for stiff ordinary differential equations,
	       Numer. Math. 33(1979), 55-68
      [RODAS3] Blom, Carmichael, Potra, Sandu, Spee, Verwer:
               Benchmarking Stiff ODE Solvers for Atmospheric Chemistry Problems II:
	       Rosenbrock Solvers,
	       Atmos. Environ. 31(1997), 3459-3472
      [ROS2]   Blom, Hundsdorfer, Spee, Verwer:
               A Second-Order Rosenbrock Method Applied to Photochemical Dispersion Problems,
	       SIAM J. Sci. Comput. 20(1999), 1456-1480
      [ROS3],
      [ROWDA3] Roche:
               Rosenbrock Methods for Differential Algebraic Equations,
               Numer. Math. 52(1988), 45-63
     */
    enum Method {
      GRK4T,     // s=4, p=4
      RODAS3,    // s=4, p=3, L-stable, stiffly accurate
      ROS2,      // s=2, p=2, L-stable
      ROS3,      // s=3, p=3, L-stable
      ROWDA3,    // s=3, p=3, L-stable
      RS,        // s=3, p=2, L-stable
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

    /*!
      set preprocessor for the right-hand side
    */
    void set_preprocessor(WMethodPreprocessRHSHelper<VECTOR>* pp) {
      preprocessor = pp;
    }
    
  protected:
    //! A=(a_{i,j})
    LowerTriangularMatrix<double> A;

    //! C=(c_{i,j})
    LowerTriangularMatrix<double> C;

    //! m=(m_i)
    Vector<double> m;

    //! coefficients for the lower order error estimator e=(e_i)
    Vector<double> e;

    //! (alpha_i)
    Vector<double> alpha_vector;

    //! (gamma_i)
    Vector<double> gamma_vector;

    //! instance of the stage equation helper
    const WMethodStageEquationHelper<VECTOR>* stage_equation_helper;

    //! instance of the RHS preprocessor
    WMethodPreprocessRHSHelper<VECTOR>* preprocessor;

    //! consistency order
    int p;
  };
}

#include <numerics/w_method.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
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
  template <class VECTOR, class IVP = AbstractIVP<VECTOR> >
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
    virtual void solve_W_stage_equation(IVP* ivp,
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
    virtual void approximate_ft(IVP* ivp,
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
				      const double tolerance)
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
  template <class VECTOR, class IVP = AbstractIVP<VECTOR> >
  class WMethod
    : public OneStepScheme<VECTOR, IVP>
  {
  public:
    /*!
      enum type for the different builtin methods

      References:
      [GRK4T]   Kaps, Rentrop:
                Generalized Runge-Kutta methods of order four with stepsize control
	        for stiff ordinary differential equations,
	        Numer. Math. 33(1979), 55-68
      [ROS3],
      [RODAS3]  Blom, Carmichael, Potra, Sandu, Spee, Verwer:
                Benchmarking Stiff ODE Solvers for Atmospheric Chemistry Problems II:
	        Rosenbrock Solvers,
	        Atmos. Environ. 31(1997), 3459-3472
      [ROS2]    Blom, Hundsdorfer, Spee, Verwer:
                A Second-Order Rosenbrock Method Applied to Photochemical Dispersion Problems,
	        SIAM J. Sci. Comput. 20(1999), 1456-1480
      [ROWDA3]  Roche:
                Rosenbrock Methods for Differential Algebraic Equations,
                Numer. Math. 52(1988), 45-63
      [ROS3P]   Lang, Verwer:
                ROS3P - An accurate third-order Rosenbrock solver designed for parabolic problems,
	        BIT 41(2001), 731--738
      [ROS3Pw]  Rang, Angermann:
                Creating new Rosenbrock methods with Maple,
                Mathematik-Bericht 2003/5, Institut fuer Mathematik, TU Clausthal
      [ROSI2P2],
      [ROSI2PW] Rang, Angermann:
                New Rosenbrock methods of order 3 for PDAEs of index 2
      [RODAS]   Hairer, Wanner:
                Solving Ordinairy Differential Equations II
		Springer Series in Computational Mathematics, Springer 1991
      [RODASP]  Steinebach:
                Order-reduction of ROW-methods for DAEs and method of lines applications
	        Technical report 1741, TH Darmstadt, 1995
     */
    enum Method {
      GRK4T,     // s=4, p=4
      RODAS3,    // s=4, p=3, index 1, L-stable, stiffly accurate
      ROS2,      // s=2, p=2, L-stable
      ROS3,      // s=3, p=3, L-stable
      ROWDA3,    // s=3, p=3, index 1, L-stable
      ROS3P,     // s=3, p=3, index 1, nonlinear PDEs
      ROS3Pw,    // s=3, p=3, index 1, PDEs
      ROSI2P2,   // s=4, p=3, index 2, L-stable, PDEs, stiffly accurate
      ROSI2PW,   // s=4, p=3, index 2, L-stable, PDEs, stiffly accurate, W-method
      RODAS,     // s=6, p=4
      RODASP     // s=6, p=4, index 1, PDEs, L-stable, stiffly accurate
    };

    /*!
      constructor from some builtin methods and a given stage equation solver
    */
    WMethod(const Method method,
	    const WMethodStageEquationHelper<VECTOR, IVP>* stage_equation_solver);
    
    /*!
      virtual destructor
    */
    virtual ~WMethod() {}

    /*!
      increment function u^{(m)} -> u^{(m+1)},
      also returns a local error estimator
    */
    void increment(IVP* ivp,
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
    
    /*!
      helper function to compute the data {A, C, m, e, alpha_vector, gamma_vector}
      from the original data {Alpha, Gamma, b, bhat},
      cf. Hairer/Wanner II.
    */
    static void transform_coefficients(const LowerTriangularMatrix<double>& Alpha,
				       const LowerTriangularMatrix<double>& Gamma,
				       const Vector<double>& b,
				       const Vector<double>& bhat,
				       LowerTriangularMatrix<double>& A,
				       LowerTriangularMatrix<double>& C,
				       Vector<double>& m,
				       Vector<double>& e,
				       Vector<double>& alpha_vector,
				       Vector<double>& gamma_vector);

    /*!
      helper function to check the algebraic order conditions of a (RO)W-method
    */
    static void check_order_conditions(const LowerTriangularMatrix<double>& Alpha,
				       const LowerTriangularMatrix<double>& Gamma,
				       const Vector<double>& b,
				       const Vector<double>& bhat,
				       const bool wmethod=false);

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


  protected:
    //! instance of the stage equation helper
    const WMethodStageEquationHelper<VECTOR, IVP>* stage_equation_helper;

    //! instance of the RHS preprocessor
    WMethodPreprocessRHSHelper<VECTOR>* preprocessor;

    //! consistency order
    int p;
  };
}

#include <numerics/w_method.cpp>

#endif

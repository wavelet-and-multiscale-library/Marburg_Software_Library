// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ITERATSOLV_H
#define _MATHTL_ITERATSOLV_H

// A collection of template-based iterative solvers for linear systems.
// 
// VECTOR: arbitrary vector class, should implement the operators +, -, * (scalar mult.)
// MATRIX: matrix class
// PREC:   preconditioner, we require apply_preconditioner(const VECTOR&, VECTOR&)

namespace FLAT
{
  // NOTE: the following routines access matrices and vectors with unsigned integer
  // indices from 0 to size()-1 !!!

  // Stationary iterative schemes (Richardson, Jacobi, Gauss-Seidel, SOR, SSOR)
  //   x^{(k+1)} = Bx^{(k)} + c

  //! Richardson iteration
  /*!
    Richardson iteration (for quadratic matrices; routine only uses APPLY())
    \param A quadratic matrix
    \param b right-hand side vector
    \param xk starting vector and solution
    \param omega relaxation parameter
    \param tol stopping tolerance (||r^{(k)}||_2 <= tol)
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  void Richardson(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
		  const double tol, const unsigned int maxiter, unsigned int& iterations);

  //! Jacobi iteration
  /*!
    Jacobi iteration (for quadratic matrices; routine only uses APPLY() and the values on the diagonal)
    \param A quadratic matrix
    \param b right-hand side vector
    \param xk starting vector and solution
    \param tol stopping tolerance (||r^{(k)}||_2 <= tol)
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  void Jacobi(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	      const double tol, const unsigned int maxiter, unsigned int& iterations);

  //! Gauss-Seidel iteration
  /*!
    Gauss-Seidel iteration (for quadratic matrices; routine does not use nonzero pattern)
    \param A quadratic matrix
    \param b right-hand side vector
    \param xk solution
    \param tol stopping tolerance (||r^{(k)}||_2 <= tol)
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  void GaussSeidel(const MATRIX &A, const VECTOR &b, VECTOR &xk,
		   const double tol, const unsigned int maxiter, unsigned int& iterations);

  //! SOR
  /*!
    SOR iteration (for quadratic matrices; routine does not use nonzero pattern)
    \param A quadratic matrix
    \param b right-hand side vector
    \param xk starting vector and solution
    \param tol stopping tolerance (||r^{(k)}||_2 <= tol)
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  void SOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const unsigned int maxiter, unsigned int& iterations);
  
  //! SSOR
  /*!
    SSOR iteration (for symmetric matrices; routine does not use nonzero pattern)
    \param A symmetric matrix
    \param b right-hand side vector
    \param xk starting vector and solution
    \param tol stopping tolerance (||r^{(k)}||_2 <= tol)
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  void SSOR(const MATRIX &A, const VECTOR &b, VECTOR &xk, const double omega,
	   const double tol, const unsigned int maxiter, unsigned int& iterations);
 
  //! conjugate gradient iteration
  /*!
    classical unpreconditioned conjugate gradient iteration
    \param A s.p.d. matrix
    \param b right-hand side vector
    \param xk starting vector and solution
    \param tol stopping tolerance
    \param iterations number of cg iterations
    \return convergence within <maxiter> iterations
  */
  template <class VECTOR, class MATRIX>
  bool CG(const MATRIX &A, const VECTOR &b, VECTOR &xk,
	  const double tol, const unsigned int maxiter, unsigned int& iterations);

  //! preconditioned conjugate gradient iteration
  /*!
    classical preconditioned conjugate gradient iteration
    \param A s.p.d. matrix
    \param b right-hand side vector
    \param P (left, s.p.d.) preconditioner
    \param xk starting vector and solution
    \param tol stopping tolerance
    \param iterations number of cg iterations
    \return convergence within <maxiter> iterations
  */
  template <class VECTOR, class MATRIX, class PREC>
  bool PCG(const MATRIX &A, const VECTOR &b, const PREC& P, VECTOR &xk,
	   const double tol, const unsigned int maxiter, unsigned int& iterations);
}

#include <numerics/iteratsolv.cpp>

#endif

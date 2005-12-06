// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_EIGENVALUES_H
#define _MATHTL_EIGENVALUES_H

// A small collection of iterative algorithms to compute
// eigenvalues of a square matrix A

// VECTOR: arbitrary vector class, should implement the operators +, -, * (scalar mult.)
// MATRIX: matrix class, we only require a routine APPLY(const VECTOR&, VECTOR&)
// PREC:   preconditioner, we again require APPLY(const VECTOR&, VECTOR&)

namespace MathTL
{
  //! von Mises power iteration
  /*!
    von Mises power iteration, cf. [Golub/van Loan], returns maximal eigenvalue in modulus
    \param A symmetric matrix
    \param xk eigenvector, starting vector (hint: choose a random vector)
    \param tol stopping tolerance
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  double PowerIteration(const MATRIX& A, VECTOR& xk,
			const double tol, const unsigned int maxit, unsigned int &iterations);
  
  //! inverse power iteration
  /*!
    Wielandt inverse power iteration, cf. Golub/van Loan,
    returns minimal eigenvalue in modulus
    \param A symmetric matrix (due to the internal CG() call)
    \param xk eigenvector, starting vector (hint: choose a random vector)
    \param tol stopping tolerance
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  double InversePowerIteration(const MATRIX& A, VECTOR& xk,
			       const double tol, const unsigned int maxit, unsigned int &iterations);

  //! inverse power iteration with spectral shift
  /*!
    Wielandt inverse power iteration (with spectral shift), cf. Golub/van Loan,
    returns minimal eigenvalue in modulus
    \param A symmetric matrix (due to the internal CG() call)
    \param xk eigenvector (1-normalized), starting vector (hint: choose a random vector)
    \param lambda shift parameter
    \param tol stopping tolerance
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  double InversePowerIteration(const MATRIX& A, VECTOR& xk, const double lambda,
			       const double tol, const unsigned int maxit, unsigned int &iterations);

  /*!
    Lanczos iteration, compute a finite subset (here the extremal points) of the spectrum of a symmetric matrix
    The iteration stops if the extremal eigenvalues have changed only below a given (relative) tolerance.
   */
  template <class MATRIX>
  void LanczosIteration(const MATRIX& A, const double tol,
			double& lambdamin, double& lambdamax, unsigned int &iterations);

  //! solve symmetric eigenvalue problem
  /*!
    Solve the symmetric eigenvalue problem
    \param A symmetric matrix
    \param evals eigenvalues
    \param evecs eigenvectors
  */
  template <class VECTOR, class MATRIX, class MATRIX2>
  void SymmEigenvalues(const MATRIX& A, VECTOR& evals, MATRIX2& evecs);

  //! spectral condition number of a symmetric matrix
  template <class MATRIX>
  double CondSymm(const MATRIX& A,
		  double tol = 1e-6, unsigned int maxit = 100);

  //! spectral condition number of a potentially nonsymmetric matrix
  template <class MATRIX>
  double CondNonSymm(const MATRIX &A,
		     double tol = 1e-6, unsigned int maxit = 100);
}

#include "numerics/eigenvalues.cpp"

#endif

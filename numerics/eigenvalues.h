// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
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
    von Mises power iteration, cf. [Golub/van Loan]
    \param A symmetric matrix
    \param xk eigenvector (1-normalized), starting vector
    \param tol stopping tolerance
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  double PowerIteration(const MATRIX& A, VECTOR& xk,
			const double tol, const unsigned int maxit, unsigned int &iterations);
  
  //! inverse power iteration
  /*!
    von Mises inverse power iteration, cf. Golub/van Loan
    \param A symmetric matrix
    \param xk eigenvector (1-normalized), starting vector
    \param tol stopping tolerance
    \param maxiter maximal number of iterations
    \param iterations iterations performed
  */
  template <class VECTOR, class MATRIX>
  double InversePowerIteration(const MATRIX& A, VECTOR& xk,
			       const double tol, const unsigned int maxit, unsigned int &iterations);

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

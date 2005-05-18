// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_DECOMPOSABLE_MATRIX_H
#define _MATHTL_DECOMPOSABLE_MATRIX_H

#include <iostream>
#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <utils/array1d.h>

// matrix norms, for convenience
#include <algebra/matrix_norms.h>

namespace MathTL
{
  /*!
    This class models finite, densely populated matrices
      M = (m_{i,j})_{0<=i<m, 0<=j<n}
    with entries from an arbitrary (scalar) class C,
    designed for the solution of linear systems
      A*x=b
    or least squares problems
      min \|A*x-b\|_2
    by the application of LU or QR factorization (m >= n).
    The factorizations is stored in place.
  */
  template <class C>
  class DecomposableMatrix
    : protected Matrix<C>
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;

    /*!
      type of factorizations
    */
    enum DecompositionType {
      none,
      LU,
      QR
    };

    /*!
      default constructor, yields zero square matrix which is empty per default
    */
    explicit DecomposableMatrix(const size_type n = 0);

    /*!
      copy constructor
    */
    DecomposableMatrix(const DecomposableMatrix<C>& M);

    /*!
      copy constructor from a Matrix
    */
    DecomposableMatrix(const Matrix<C>& M);

    /*!
      row dimension
    */
    const size_type row_dimension() const;

    /*!
      column dimension
    */
    const size_type column_dimension() const;

    /*!
      size as an STL-compatible container for matrix entries
    */
    const size_type size() const;

    /*!
      perform (or revert) decomposition
    */
    void decompose(DecompositionType d = LU);

    /*!
      solve the linear system A*x=b via a given decomposition
      (x will be scaled appropriately)
    */
    void solve(const Vector<C>& b, Vector<C>& x,
	       DecompositionType d = LU);

     /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 10,
	       const unsigned int precision = 3) const;

 protected:
    // decomposition type
    DecompositionType decomposition;

    // storage for row equilibration matrix entries (-> LU factorization)
    Vector<double> D;

    // storage for permuation matrix (-> LU factorization)
    Array1D<size_type> P;

    // LU decomposition (row equilibration, row pivoting, in place (up to D and P))
    void LU_decomposition();

    // QR decomposition (in place)
    void QR_decomposition();

    // revert decomposition
    void revert_decomposition();
  };

  /*!
    Matlab-style stream output for decomposable matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const DecomposableMatrix<C>& M);
}

// include implementation of inline functions
#include <numerics/decomposable_matrix.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MATRIX_H
#define _MATHTL_MATRIX_H

#include <iostream>
#include <algebra/vector.h>
#include <algebra/symmetric_matrix.h>

// matrix norms, for convenience
#include <algebra/matrix_norms.h>

namespace MathTL
{
  /*!
    This class models finite, densely populated matrices
      M = (m_{i,j})_{0<=i<m, 0<=j<n}
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
  */
  template <class C>
  class Matrix
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;
    
    /*!
      default constructor, yields zero square matrix which is empty per default
    */
    explicit Matrix(const size_type n = 0);

    /*!
      copy constructor
    */
    Matrix(const Matrix<C>& M);

    /*!
      copy constructor from symmetric matrices
    */
    Matrix(const SymmetricMatrix<C>& M);
    
    /*!
      construct m*n rectangular matrix
    */
    Matrix(const size_type row_dimension, const size_type column_dimension);

    /*!
      Construct matrix from a string holding its entries, separated
      by a blank.
      \param row_dimension row dimension
      \param column_dimension column dimension
      \param str input string
      \param byrow indicates whether coefficients are stored row by row in the stream
     */
    Matrix(const size_type row_dimension,
	   const size_type column_dimension,
	   const char* str,
	   const bool byrow = true);

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
      resize matrix and initialize with zero
    */
    void resize(const size_type rows, const size_type columns);

    /*!
      (estimate for the) memory consumption in bytes
    */
    const size_type memory_consumption() const;

    /*!
      return true if matrix is empty (cf. STL containers)
    */
    bool empty() const;

    /*!
      read-only access to a matrix entry
    */
    const C operator () (const size_type row, const size_type column) const;

    /*!
      read-write access to a matrix entry
    */
    C& operator () (const size_type row, const size_type column);

    /*!
      equality test with another matrix
    */
    template <class C2>
    bool operator == (const Matrix<C2>& M) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const Matrix<C2>& M) const;

    /*!
      assignment from another matrix
    */
    Matrix<C>& operator = (const Matrix<C>& M);

    /*!
      swap entries of two matrices
    */
    void swap (Matrix<C>& M);
    
    /*!
      matrix-vector multiplication Mx = (*this) * x;
      we assume that the vector Mx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& Mx) const;

    /*!
      transposed matrix-vector multiplication Mtx = (*this)^T * x;
      we assume that the vector Mtx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& Mtx) const;

    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 5,
	       const unsigned int precision = 2) const;

  protected:
    /*!
      internal storage of densely populated matrices is just an
      appropriately sized vector, which holds the matrix entries
      in row major ordering
    */
    Vector<C> entries_;

    /*!
      row dimension
    */
    size_type rowdim_;

    /*!
      column dimension
    */
    size_type coldim_;
  };

  /*!
    swap the entries of two matrices
  */
  template <class C>
  void swap(Matrix<C>& M1, Matrix<C>& M2);

  /*!
    Matlab-style stream output for dense matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const Matrix<C>& M);
}

// include implementation of inline functions
#include <algebra/matrix.cpp>

#endif

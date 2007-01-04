// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_TRIDIAGONAL_MATRIX_H
#define _MATHTL_TRIDIAGONAL_MATRIX_H

#include <iostream>
#include <algebra/vector.h>

namespace MathTL
{
  /*!
    This class models tridiagonal matrices
      M = (m_{i,j})_{0<=i,j<n}
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
    
    Internally, the entries are stored as three vectors.
  */
  template <class C>
  class TridiagonalMatrix
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;

    /*!
      default constructor, yields zero tridiagonal matrix which is 1x1 by default
    */
    explicit TridiagonalMatrix(const size_type n = 1);

    /*!
      copy constructor
    */
    TridiagonalMatrix(const TridiagonalMatrix<C>& M);

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
    void resize(const size_type rows);

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
      read-only access to a matrix entry
    */
    const C get_entry(const size_type row, const size_type column) const;

    /*!
      read-write access to a matrix entry
    */
    C& operator () (const size_type row, const size_type column);
    
    /*!
      write access to a matrix entry
    */
    void set_entry(const size_type row, const size_type column, const C value);

    /*!
      equality test with another matrix
    */
    template <class C2>
    bool operator == (const TridiagonalMatrix<C2>& M) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const TridiagonalMatrix<C2>& M) const;

    /*!
      assignment from another tridiagonal matrix
    */
    TridiagonalMatrix<C>& operator = (const TridiagonalMatrix<C>& M);

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
      a: lower diagonal
      b: main diagonal
      c. upper diagonal
    */
    Vector<C> a_, b_, c_;

    /*!
      row dimension == column dimension
    */
    size_type rowdim_;
  };

  /*!
    Matlab-style stream output for symmetric matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const TridiagonalMatrix<C>& M);
}

// include implementation of inline functions
#include <algebra/tridiagonal_matrix.cpp>

#endif

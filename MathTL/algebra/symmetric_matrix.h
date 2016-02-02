// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SYMMETRIC_MATRIX_H
#define _MATHTL_SYMMETRIC_MATRIX_H

#include <iostream>
#include <algebra/triangular_matrix.h>

namespace MathTL
{
  /*!
    This class models finite, densely populated, symmetric matrices
      M = (m_{i,j})_{0<=i<m, 0<=j<n}
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
    
    Internally, we only store the lower triangular part.
  */
  template <class C>
  class SymmetricMatrix
    : public LowerTriangularMatrix<C>
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename LowerTriangularMatrix<C>::size_type size_type;

    /*!
      default constructor, yields zero square matrix which is empty per default
    */
    explicit SymmetricMatrix(const size_type n = 0);

    /*!
      copy constructor
    */
    SymmetricMatrix(const SymmetricMatrix<C>& M);

    /*!
      Construct n*n matrix from a string holding its entries, separated
      by a blank. Please note that we internally only store the lower triangular
      part of the matrix, so the string only has to contain those entries.
      \param n row and column dimension
      \param str input string
      \param byrow indicates whether coefficients are stored row by row in the stream
     */
    SymmetricMatrix(const size_type n,
		    const char* str,
		    const bool byrow = true);

    /*!
      resize matrix and initialize with zero
    */
    void resize(const size_type n);

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
      (you should not use this, since it destroys the matrix structure!)
    */
    C& operator () (const size_type row, const size_type column);

    /*!
      write access to a matrix entry
      (you should not use this, since it destroys the matrix structure!)
    */
    void set_entry(const size_type row, const size_type column, const C value);

    /*!
      equality test with another matrix
    */
    template <class C2>
    bool operator == (const SymmetricMatrix<C2>& M) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const SymmetricMatrix<C2>& M) const;

    /*!
      assignment from another symmetric matrix
    */
    SymmetricMatrix<C>& operator = (const SymmetricMatrix<C>& M);

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
      (this is identical to apply(), since we have a symmetric matrix)
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
  };

  /*!
    Matlab-style stream output for symmetric matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const SymmetricMatrix<C>& M);
}

// include implementation of inline functions
#include <algebra/symmetric_matrix.cpp>

#endif

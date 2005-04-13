// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SPARSE_MATRIX_H
#define _MATHTL_SPARSE_MATRIX_H

#include <iostream>

#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/symmetric_matrix.h>
#include <algebra/triangular_matrix.h>

// matrix norms, for convenience
#include <algebra/matrix_norms.h>

namespace MathTL
{
  /*!
    This class models finite, sparsely populated matrices
      M = (m_{i,j})_{0<=i<m, 0<=j<n}
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
    The internal representation is CRS (compressed row storage), see [N].

    Reference:
    [N] http://www.netlib.org/linalg/html_templates/node91.html
  */
  template <class C>
  class SparseMatrix
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;

    /*!
      default constructor, yields zero square matrix which is empty per default
    */
    explicit SparseMatrix(const size_type n = 0);

    /*!
      destructor
    */
    ~SparseMatrix();
    
    /*!
      row dimension
    */
    const size_type row_dimension() const;

    /*!
      column dimension
    */
    const size_type column_dimension() const;

    /*!
      size as an STL-compatible container for matrix entries,
      number of nonzero entries
    */
    const size_type size() const;

    /*!
      resize matrix and initialize with zero
    */
    void resize(const size_type rows, const size_type columns);

    /*!
      read-only access to a matrix entry
    */
    const C get_entry(const size_type row, const size_type column) const;

    /*!
      read-write access to a matrix entry
    */
    void set_entry(const size_type row, const size_type column, const C value);

    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 5,
	       const unsigned int precision = 2) const;

  protected:
    /*!
      storage for the matrix entries,
      entries_[r] contains the elements in row r
    */
    C** entries_;

    /*!
      storage for the indices of nontrivial entries,
      indices_[r] contains the nonzero indices in row r, where
      * indices_[r][0] is the number of indices k
      * indices_[r][1],...,indices_[r][k] are the indices themselves
    */
    size_type** indices_;

    /*!
      row dimension
    */
    size_type rowdim_;

    /*!
      column dimension
    */
    size_type coldim_;

    /*!
      release memory
    */
    void clear();
  };

  /*!
    Matlab-style stream output as a dense matrix
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const SparseMatrix<C>& M);
}

// include implementation of inline functions
#include <algebra/sparse_matrix.cpp>

#endif

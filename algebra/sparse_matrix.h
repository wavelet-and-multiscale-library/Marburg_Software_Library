// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SPARSE_MATRIX_H
#define _MATHTL_SPARSE_MATRIX_H

#include <iostream>
#include <list>

#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/matrix_block.h>
#include <algebra/symmetric_matrix.h>
#include <algebra/triangular_matrix.h>
#include <algebra/infinite_vector.h>

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
    : public MatrixBlock<C>
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;

    /*!
      default constructor, yields zero square matrix
    */
    explicit SparseMatrix(const size_type n = 1);

    /*!
      copy constructor
    */
    SparseMatrix(const SparseMatrix<C>& M);

    /*!
      construct m*n rectangular zero matrix
    */
    SparseMatrix(const size_type row_dimension, const size_type column_dimension);

    /*!
      destructor
    */
    ~SparseMatrix();
    
    //! clone the sparse matrix (requirement from MatrixBlock)
    MatrixBlock<C>* clone() const;

    //! transpose the sparse matrix (requirement from MatrixBlock)
    MatrixBlock<C>* clone_transposed() const;

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
      return true if matrix is empty (cf. STL containers)
    */
    bool empty() const;

    /*!
      resize matrix and initialize with zero
    */
    void resize(const size_type rows, const size_type columns);

    /*!
      number of nonzero entries a given row
    */
    const size_type entries_in_row(const size_type row) const;
    
    /*!
      read access to the column index of the n-th nontrivial element in a given row
    */
    const size_type get_nth_index(const size_type row, const size_type n) const;

    /*!
      read access to the n-th nontrivial element in a given row
    */
    const C get_nth_entry(const size_type row, const size_type n) const;

    /*!
      read-only access to a single matrix entry
    */
    const C get_entry(const size_type row, const size_type column) const;

    /*!
      read access to a subblock
    */
    template <class MATRIX>
    void get_block(const size_type firstrow, const size_type firstcolumn,
		   const size_type rows, const size_type columns,
		   MATRIX& M) const;
    
    /*!
      read access to an entire row;
      offset leads to a shift of v (with respect to the original column indices)
    */
    void get_row(const size_type row,
		 InfiniteVector<C, size_type>& v,
		 const size_type offset = 0) const;

    /*!
      write access to a matrix entry
    */
    void set_entry(const size_type row, const size_type column, const C value);

    /*!
      write access to a complete row
      (this std::list version is especially suitable for situations where the nonzero
      pattern is not a priorily known)
    */
    void set_row(const size_type row,
		 const std::list<size_type>& indices,
		 const std::list<C>& entries);

    /*!
      write access to a subblock;
      if the mirror flag is set, rows and columns of the block are mirrored before writing
    */
    template <class MATRIX>
    void set_block(const size_type firstrow, const size_type firstcolumn,
		   const MATRIX& M,
		   const bool mirror = false);

    /*!
      assignment from another sparse matrix
    */
    SparseMatrix<C>& operator = (const SparseMatrix<C>& M);

    /*!
      in place scaling *this *= s
    */
    void scale(const C s);

    /*!
      yields an n-by-n diagonal matrix
    */
    void diagonal(const size_type n, const C diag);

    /*!
      matrix-vector multiplication Mx = (*this) * x;
      we assume that the vector Mx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& Mx) const;

    //! special version for Vector<C> (requirement from MatrixBlock)
    void apply(const Vector<C>& x, Vector<C>& Mx) const;

    /*!
      transposed matrix-vector multiplication Mtx = (*this)^T * x;
      we assume that the vector Mtx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& Mtx) const;

    //! special version for Vector<C> (requirement from MatrixBlock)
    void apply_transposed(const Vector<C>& x, Vector<C>& Mtx) const;
    
    /*!
      set all values with modulus below a threshold to zero
      (fabs<C> should exist)
    */
    void compress(const double eta = 1e-16);

    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 8,
	       const unsigned int precision = 3) const;

    /*!
      write matlab sparse matrix format.
    */
    void matlab_output (const char *file, char *Matrixname, int binary) const;

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
      deallocate all memory
    */
    void kill();
  };

  /*!
    matrix-matrix difference M-N
  */
  template <class C>
  SparseMatrix<C> operator - (const SparseMatrix<C>& M, const SparseMatrix<C>& N);

  /*!
    matrix-matrix multiplication M*N
  */
  template <class C>
  SparseMatrix<C> operator * (const SparseMatrix<C>& M, const SparseMatrix<C>& N);

  /*!
    transpose of a matrix
  */
  template <class C>
  SparseMatrix<C> transpose(const SparseMatrix<C>& M);

  /*!
    Matlab-style stream output as a dense matrix
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const SparseMatrix<C>& M);
}

// include implementation of inline functions
#include <algebra/sparse_matrix.cpp>

#endif

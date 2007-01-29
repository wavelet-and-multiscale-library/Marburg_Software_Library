// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_FIXED_MATRIX_H
#define _MATHTL_FIXED_MATRIX_H

#include <iostream>
#include <algebra/fixed_vector.h>

// matrix norms, for convenience
#include <algebra/matrix_norms.h>

namespace MathTL
{
  /*!
    This class models finite, densely populated matrices
      M = (m_{i,j})_{0<=i<m, 0<=j<n}
    with entries from an arbitrary (scalar) class C, where the matrix size
    is a priori known.
  */
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM = ROW_DIM>
  class FixedMatrix
  {
    public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename FixedVector<C,ROW_DIM*COL_DIM>::size_type size_type;
    
    /* constructors and destructors ****************************************/
    /*!
      default constructor, yields zero square matrix which is empty per default
    */
    explicit FixedMatrix();

    /*!
      copy constructor
    */
    FixedMatrix(const FixedMatrix<C, ROW_DIM, COL_DIM>& M);

    /*!
      Construct matrix from a string holding its entries, separated by a blank.
      \param str input string
      \param byrow indicates whether coefficients are stored row by row in the stream
     */
    FixedMatrix(const char* str,
                const bool byrow = true);

    /*!
      release allocated memory
    */
    ~FixedMatrix();

    /* member functions and operators **************************************/
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
    const C get_entry(const size_type row, const size_type column) const;

    /*!
      read-only access to a matrix entry
    */
    const C operator () (const size_type row, const size_type column) const;

    /*!
      write access to a matrix entry
    */
    void set_entry(const size_type row, const size_type column, const C value);

    /*!
      read-write access to a matrix entry
    */
    C& operator () (const size_type row, const size_type column);
    
    /*!
      equality test with another fixed matrix
    */
    template <class C2>
    bool operator == (const FixedMatrix<C2, ROW_DIM, COL_DIM>& M) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const FixedMatrix<C2, ROW_DIM, COL_DIM>& M) const;

    /*!
      assignment from another fixed matrix
    */
    FixedMatrix<C, ROW_DIM, COL_DIM>& operator = (const FixedMatrix<C, ROW_DIM, COL_DIM>& M);

    /*!
      in place scaling *this *= s
    */
    void scale(const C s);

    /*!
      in place scaling *this *= s
    */
    void operator *= (const C s);

    /*!
      in place addition *this += M
    */
    void add(const FixedMatrix<C, ROW_DIM, COL_DIM>& M);

    /*!
      in place addition *this += M
    */
    void operator += (const FixedMatrix<C, ROW_DIM, COL_DIM>& M);

    /*!
      in place subtraction *this -= M
    */
    void operator -= (const FixedMatrix<C, ROW_DIM, COL_DIM>& M);

    /*!
      matrix-vector multiplication Mx = (*this) * x;
      we assume that the vector Mx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& Mx) const;

    //! special version for FixedVector<C, DIM> (requirement from MatrixBlock)
    void apply(const FixedVector<C, COL_DIM>& x, FixedVector<C, ROW_DIM>& Mx) const;

    /*!
      transposed matrix-vector multiplication Mtx = (*this)^T * x;
      we assume that the vector Mtx has the correct size and
      is not identical to x
    */
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& Mtx) const;

    //! special version for FixedVector<C, DIM> (requirement from MatrixBlock)
    void apply_transposed(const FixedVector<C, ROW_DIM>& x, FixedVector<C, COL_DIM>& Mtx) const;

    /*!
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
               const unsigned int tabwidth = 10,
               const unsigned int precision = 3) const;


    /* private data fields *************************************************/
    private:
    /*!
      internal storage of densely populated matrices is just an
      appropriately sized vector, which holds the matrix entries
      in row major ordering
    */
    FixedVector<C, ROW_DIM*COL_DIM> entries_;
  };


  /*!
    matrix-matrix difference M-N
  */
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  FixedMatrix<C, ROW_DIM, COL_DIM> operator - (const FixedMatrix<C, ROW_DIM, COL_DIM>& M, const FixedMatrix<C, ROW_DIM, COL_DIM>& N);

  /*!
    matrix-matrix multiplication M*N
  */
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM, unsigned int COL_DIM2>
  FixedMatrix<C, ROW_DIM, COL_DIM> operator * (const FixedMatrix<C, ROW_DIM, COL_DIM>& M, const FixedMatrix<C, COL_DIM, COL_DIM2>& N);

  /*!
    transpose of a matrix
  */
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  FixedMatrix<C, COL_DIM, ROW_DIM> transpose(const FixedMatrix<C, ROW_DIM, COL_DIM>& M);

  /*!
    Matlab-style stream output for dense matrices
  */
  template <class C, unsigned int ROW_DIM, unsigned int COL_DIM>
  std::ostream& operator << (std::ostream& os, const FixedMatrix<C, ROW_DIM, COL_DIM>& M);

}

// include implementation of inline functions
#include <algebra/fixed_matrix.cpp>

#endif // _MATHTL_FIXED_MATRIX_H

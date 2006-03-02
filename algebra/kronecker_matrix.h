// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_KRONECKER_MATRIX_H
#define _MATHTL_KRONECKER_MATRIX_H

#include <iostream>
#include <algebra/vector.h>

// matrix norms, for convenience
#include <algebra/matrix_norms.h>

namespace MathTL
{
  /*!
    This class models finite matrices stemming from a Kronecker product
      kron(A,B) = (a_{1,1}B ... a_{1,n}B)
                  (.                    )
                  (a_{m,1}B ... a_{m,n}B)
    of two (arbitrary) matrices A and B.
    The entries of M stem from an arbitrary (scalar) class C.
  */
  template <class C, class MATRIX1, class MATRIX2>
  class KroneckerMatrix
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;
    
    /*!
      default constructor from A and B
    */
    explicit KroneckerMatrix(const MATRIX1& A, const MATRIX2& B);

//   //! given a matrix A, provide A^TA as a matrix
//   template <class MATRIX>
//   class AtrA
//   {
//   public:
//     /*!
//       size type
//      */
//     typedef typename MATRIX::size_type size_type;

//     //! default constructor
//     AtrA(const MATRIX& A) : A_(A) {}

//     //! row dimension
//     const size_type row_dimension() const { return A_.column_dimension(); }
    
//     //! column dimension
//     const size_type column_dimension() const { return A_.column_dimension(); }
    
//     //! apply A^TA
//     template <class VECTOR>
//     void apply(const VECTOR& x, VECTOR& AtrAx) const;
    
//     //! apply (A^TA)^T = A^TA
//     template <class VECTOR>
//     void apply_transposed(const VECTOR& x, VECTOR& AtrAx) const { apply(x, AtrAx); }

//   private:
//     const MATRIX& A_;
//   };


//     /*!
//       copy constructor
//     */
//     Matrix(const Matrix<C>& M);


//     /*!
//       construct m*n rectangular matrix
//     */
//     Matrix(const size_type row_dimension, const size_type column_dimension);

//     /*!
//       Construct matrix from a string holding its entries, separated
//       by a blank.
//       \param row_dimension row dimension
//       \param column_dimension column dimension
//       \param str input string
//       \param byrow indicates whether coefficients are stored row by row in the stream
//      */
//     Matrix(const size_type row_dimension,
// 	   const size_type column_dimension,
// 	   const char* str,
// 	   const bool byrow = true);

//     /*!
//       row dimension
//     */
//     const size_type row_dimension() const;

//     /*!
//       column dimension
//     */
//     const size_type column_dimension() const;

//     /*!
//       size as an STL-compatible container for matrix entries
//     */
//     const size_type size() const;
    
//     /*!
//       resize matrix and initialize with zero
//     */
//     void resize(const size_type rows, const size_type columns);

//     /*!
//       resize the matrix, prepare an n-by-n diagonal matrix
//     */
//     void diagonal(const size_type n, const C diag);

//     /*!
//       (estimate for the) memory consumption in bytes
//     */
//     const size_type memory_consumption() const;

//     /*!
//       return true if matrix is empty (cf. STL containers)
//     */
//     bool empty() const;

//     /*!
//       read-only access to a matrix entry
//     */
//     const C operator () (const size_type row, const size_type column) const;

//     /*!
//       read-only access to a matrix entry
//     */
//     const C get_entry(const size_type row, const size_type column) const;

//     /*!
//       read-write access to a matrix entry
//     */
//     C& operator () (const size_type row, const size_type column);
    
//     /*!
//       write access to a matrix entry
//     */
//     void set_entry(const size_type row, const size_type column, const C value);

//     /*!
//       write access to a subblock;
//       if the mirror flag is set, rows and columns of the block are mirrored before writing
//     */
//     template <class MATRIX>
//     void set_block(const size_type firstrow, const size_type firstcolumn,
// 		   const MATRIX& M,
// 		   const bool mirror = false);

//     /*!
//       equality test with another matrix
//     */
//     template <class C2>
//     bool operator == (const Matrix<C2>& M) const;

//     /*!
//       non-equality test
//     */
//     template <class C2>
//     bool operator != (const Matrix<C2>& M) const;

//     /*!
//       assignment from another matrix
//     */
//     Matrix<C>& operator = (const Matrix<C>& M);

//     /*!
//       swap entries of two matrices
//     */
//     void swap (Matrix<C>& M);

//     /*!
//       create mirrored matrix (flip row and column numbers)
//     */
//     void mirror (Matrix<C>& M) const;
    
//     /*!
//       in place scaling *this *= s
//     */
//     void scale(const C s);

//     /*!
//       matrix-vector multiplication Mx = (*this) * x;
//       we assume that the vector Mx has the correct size and
//       is not identical to x
//     */
//     template <class VECTOR>
//     void apply(const VECTOR& x, VECTOR& Mx) const;

//     /*!
//       transposed matrix-vector multiplication Mtx = (*this)^T * x;
//       we assume that the vector Mtx has the correct size and
//       is not identical to x
//     */
//     template <class VECTOR>
//     void apply_transposed(const VECTOR& x, VECTOR& Mtx) const;

//     /*!
//       set all values with modulus below a threshold to zero
//       (fabs<C> should exist)
//     */
//     void compress(const double eta = 1e-16);

//     /*!
//       stream output with user-defined tabwidth and precision
//       (cf. deal.II)
//     */
//     void print(std::ostream& os,
// 	       const unsigned int tabwidth = 10,
// 	       const unsigned int precision = 3) const;

  protected:
    MATRIX1 A;
    MATRIX2 B;
  };

  /*!
    Matlab-style stream output for Kronecker matrices
  */
  template <class C, class MATRIX1, class MATRIX2>
  std::ostream& operator << (std::ostream& os, const KroneckerMatrix<C,MATRIX1,MATRIX2>& M);
}

// include implementation of inline functions
#include <algebra/kronecker_matrix.cpp>

#endif

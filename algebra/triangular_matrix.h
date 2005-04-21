// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_TRIANGULAR_MATRIX_H
#define _MATHTL_TRIANGULAR_MATRIX_H

#include <iostream>
#include <algebra/vector.h>

// matrix norms, for convenience
#include <algebra/matrix_norms.h>

namespace MathTL
{
  /*!
    This class models finite, densely populated, lower triangular matrices
      M = (m_{i,j})_{0<=i<m, 0<=j<n}
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
    The class also works for nonquadratic matrices, where the triangular
    shape is appropriately clipped or filled up with zeros.
  */
  template <class C>
  class LowerTriangularMatrix
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;

    /*!
      default constructor, yields zero square matrix which is empty per default
    */
    explicit LowerTriangularMatrix(const size_type n = 0);

    /*!
      copy constructor
    */
    LowerTriangularMatrix(const LowerTriangularMatrix<C>& M);

    /*!
      construct m*n rectangular matrix
    */
    LowerTriangularMatrix(const size_type rows, const size_type columns);

    /*!
      Construct matrix from a string holding its entries, separated
      by a blank.
      \param row_dimension row dimension
      \param column_dimension column dimension
      \param str input string
      \param byrow indicates whether coefficients are stored row by row in the stream
     */
    LowerTriangularMatrix(const size_type rows,
			  const size_type columns,
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
    bool operator == (const LowerTriangularMatrix<C2>& M) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const LowerTriangularMatrix<C2>& M) const;

    /*!
      assignment from another lower triangular matrix
    */
    LowerTriangularMatrix<C>& operator = (const LowerTriangularMatrix<C>& M);

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
      return number of nonzero entries in a triangular matrix
    */
    static size_type triangle_size(const size_type rows, const size_type columns);

    /*!
      return index of an entry in a lower triangular (row major storage) or upper
      triangular (column major storage) or symmetric matrix
    */
    static size_type triangle_index(const size_type row,
				    const size_type column,
				    const size_type rowdim,
				    const size_type coldim);
   
    /*!
      internal storage of lower triangular matrices is just an
      appropriately sized vector, which holds the matrix entries
      in row major ordering with decreasing width
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
    Matlab-style stream output for lower triangular matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const LowerTriangularMatrix<C>& M);


  /*!
    This class models finite, densely populated, upper triangular matrices
      M = (m_{i,j})_{0<=i<m, 0<=j<n}
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
    The class also works for nonquadratic matrices, where the triangular
    shape is appropriately clipped or filled up with zeros.
  */
  template <class C>
  class UpperTriangularMatrix
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
    explicit UpperTriangularMatrix(const size_type n = 0);

    /*!
      copy constructor
    */
    UpperTriangularMatrix(const UpperTriangularMatrix<C>& M);

    /*!
      construct m*n rectangular matrix
    */
    UpperTriangularMatrix(const size_type rows, const size_type columns);

    /*!
      Construct matrix from a string holding its entries, separated
      by a blank.
      \param row_dimension row dimension
      \param column_dimension column dimension
      \param str input string
      \param byrow indicates whether coefficients are stored row by row in the stream
     */
    UpperTriangularMatrix(const size_type rows,
			  const size_type columns,
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
      resize matrix and initialize with zero
    */
    void resize(const size_type rows, const size_type columns);

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
    bool operator == (const UpperTriangularMatrix<C2>& M) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const UpperTriangularMatrix<C2>& M) const;

    /*!
      assignment from another upper triangular matrix
    */
    UpperTriangularMatrix<C>& operator = (const UpperTriangularMatrix<C>& M);

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
  };

  /*!
    Matlab-style stream output for upper triangular matrices
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const UpperTriangularMatrix<C>& M);
}

// include implementation of inline functions
#include <algebra/triangular_matrix.cpp>

#endif

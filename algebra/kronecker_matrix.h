// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_KRONECKER_MATRIX_H
#define _MATHTL_KRONECKER_MATRIX_H

#include <iostream>
#include <algebra/vector.h>
#include <algebra/matrix_block.h>
#include <algebra/matrix.h>

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
    : public MatrixBlock<C>
  {
  public:
    //! type of indexes and size type (cf. STL containers)
    typedef typename Vector<C>::size_type size_type;
    
    //! default constructor from A and B
    explicit KroneckerMatrix(const MATRIX1& A, const MATRIX2& B,
			     const double factor = 1.0);
    
    //! copy constructor
    KroneckerMatrix(const KroneckerMatrix<C,MATRIX1,MATRIX2>& M);

    //! clone the Kronecker matrix (requirement from MatrixBlock)
    MatrixBlock<C>* clone() const;

    //! transpose the Kronecker matrix (requirement from MatrixBlock)
    MatrixBlock<C>* clone_transposed() const;

    //! row dimension
    const size_type row_dimension() const;

    //! column dimension
    const size_type column_dimension() const;

    //! return true if matrix is empty (cf. STL containers)
    bool empty() const;

    //! read-only access to a matrix entry
    const C operator () (const size_type row, const size_type column) const;

    //! read-only access to a matrix entry
    const C get_entry(const size_type row, const size_type column) const;

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
      stream output with user-defined tabwidth and precision
      (cf. deal.II)
    */
    void print(std::ostream& os,
	       const unsigned int tabwidth = 10,
	       const unsigned int precision = 3) const;

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

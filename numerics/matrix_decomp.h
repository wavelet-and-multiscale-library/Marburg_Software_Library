// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MATRIX_DECOMP_H
#define _MATHTL_MATRIX_DECOMP_H

#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/triangular_matrix.h>
#include <algebra/symmetric_matrix.h>

namespace MathTL
{
  //! Cholesky decomposition
  /*!
    Cholesky decomposition, construct L such that
      A = LL^T
    \param A s.p.d. matrix with real entries
    \param L result matrix
    \return successful decomposition (i.e. A is really s.p.d.)
  */
  template <class C>
  bool CholeskyDecomposition(const SymmetricMatrix<C>& A,
			     LowerTriangularMatrix<C>& L);

  /*!
    QU decomposition, construct unitary Q and "upper triangular" U such that
      A = QU

    References:
      JAMA
      Stoer, Numerische Mathematik I
    \param A arbitrary matrix
  */
  template <class C>
  class QUDecomposition
  {
  public:
    //! default constructor, perform QU decomposition
    QUDecomposition(const Matrix<C>& A);

    //! determine whether A has full rank
    bool hasFullRank() const;

    //! return U
    void getU(UpperTriangularMatrix<C>& U) const;

    //! return Q
    void getQ(Matrix<C>& Q) const;

    /*!
      After the QU decomposition, solve the linear system Ax = b,
      where we assume that b is in the range of A.
      The vector x will be resized properly
    */
    void solve(const Vector<C>& b, Vector<C>& x) const;

    /*!
      After the QU decomposition, compute the inverse of A.
    */
    void inverse(Matrix<C>& AInv) const;

  protected:
    //! m = A.row_dimension(), n = A.column_dimension()
    typename Matrix<C>::size_type rowdim_, coldim_;
    //! storage for the decomposition
    Matrix<C> QU_;
    //! storage for diag(U)
    Vector<C> Udiag_;
  };

  //! singular value decomposition
  /*!
    Compute the singular value decomposition of an m-by-n matrix A, in the form
      A = U*S*V
    where U is an m-by-n orthogonal matrix, S is an n-by-n diagonal matrix
    and V is an n-by-n orthogonal matrix.
    The singular values sigma_k = S(k,k) are ordered such that
      sigma_0 >= sigma_1 >= ... >= sigma_{n-1}
  */
  template <class C>
  class SVD
  {
  public:
    //! default constructor, perform SVD
    SVD(const Matrix<C>& A);

    //! return U
    void getU(Matrix<C>& U) const;

    //! return U*S
    void getUS(Matrix<C>& US) const;

    //! return singular values diag(S)
    void getS(Vector<C>& S) const { S = Sdiag_; }

    //! return V
    void getV(Matrix<C>& V) const;

  protected:
    //! m=A.rowdim(), n=A.coldim()
    typename Matrix<C>::size_type rowdim_, coldim_;
    //! storage for the decomposition
    Matrix<C> U_, V_;
    //! storage for diag(S)
    Vector<C> Sdiag_;
  };

  /*!
    Convert a given matrix A into Hessenberg form B=U^T*A*U, where U is orthogonal.
    We do not store the U-matrix, since this routine is the main ingredient
    in the numerical solution of the nonsymmetric eigenvalue problem.
  */
  template <class C>
  class Hessenberg
    : public Matrix<C>
  {
  public:
    //! default constructor, perform the Hessenberg transformation
    Hessenberg(const Matrix<C>& A);
  };
}

#include <numerics/matrix_decomp.cpp>

#endif

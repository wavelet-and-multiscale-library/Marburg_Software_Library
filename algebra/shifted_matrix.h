// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SHIFTED_MATRIX_H
#define _MATHTL_SHIFTED_MATRIX_H

#include <cassert>

namespace MathTL
{
  //! given a square matrix A, provide A-lambda*I as a matrix
  template <class MATRIX>
  class ShiftedMatrix
  {
  public:
    /*!
      size type
     */
    typedef typename MATRIX::size_type size_type;

    //! (default) constructor
    ShiftedMatrix(const MATRIX& A, const double lambda = 0)
      : A_(A), lambda_(lambda)
    {
      assert(A.row_dimension() == A.column_dimension());
    }

    //! set spectral parameter
    void set_lambda(const double lambda) { lambda_ = lambda; }

    //! row dimension
    const size_type row_dimension() const { return A_.column_dimension(); }
    
    //! column dimension
    const size_type column_dimension() const { return A_.column_dimension(); }
    
    //! apply A-lambda*I
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& result) const;
    
    //! apply (A-lambda*I)^T = A^T-lambda*I
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& result) const;
    
  private:
    const MATRIX& A_;
    double lambda_;
  };

  template <class MATRIX> template <class VECTOR>
  void ShiftedMatrix<MATRIX>::apply(const VECTOR& x, VECTOR& result) const
  {
    A_.apply(x, result);
    result.add(-lambda_, x);
  }

  template <class MATRIX> template <class VECTOR>
  void ShiftedMatrix<MATRIX>::apply_transposed(const VECTOR& x, VECTOR& result) const
  {
    A_.apply_transposed(x, result);
    result.add(-lambda_, x);
  }
}

#endif

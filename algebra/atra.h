// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ATRA_H
#define _MATHTL_ATRA_H

namespace MathTL
{
  //! given a matrix A, provide A^TA as a matrix
  template <class MATRIX>
  class AtrA
  {
  public:
    /*!
      size type
     */
    typedef typename MATRIX::size_type size_type;

    //! default constructor
    AtrA(const MATRIX& A) : A_(A) {}

    //! row dimension
    const size_type row_dimension() const { return A_.column_dimension(); }
    
    //! column dimension
    const size_type column_dimension() const { return A_.column_dimension(); }
    
    //! apply A^TA
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& AtrAx) const;
    
    //! apply (A^TA)^T = A^TA
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& AtrAx) const { apply(x, AtrAx); }

  private:
    const MATRIX& A_;
  };
  
  template <class MATRIX> template <class VECTOR>
  void AtrA<MATRIX>::apply(const VECTOR& x, VECTOR& AtrAx) const
  {
    VECTOR y(A_.row_dimension(), false);
    A_.apply(x, y);
    A_.apply_transposed(y, AtrAx);
  }

  //! given a matrix A, provide AA^T as a matrix
  template <class MATRIX>
  class AAtr
  {
  public:
    /*!
      size type
     */
    typedef typename MATRIX::size_type size_type;

    //! default constructor
    AAtr(const MATRIX& A) : A_(A) {}

    //! row dimension
    const size_type row_dimension() const { return A_.row_dimension(); }
    
    //! column dimension
    const size_type column_dimension() const { return A_.row_dimension(); }
    
    //! apply AA^T
    template <class VECTOR>
    void apply(const VECTOR& x, VECTOR& AAtrx) const;
    
    //! apply (AA^T)^T = AA^T
    template <class VECTOR>
    void apply_transposed(const VECTOR& x, VECTOR& AAtrx) const { apply(x, AAtrx); }

  private:
    const MATRIX& A_;
  };
  
  template <class MATRIX> template <class VECTOR>
  void AAtr<MATRIX>::apply(const VECTOR& x, VECTOR& AAtrx) const
  {
    VECTOR y(A_.column_dimension(), false);
    A_.apply_transposed(x, y);
    A_.apply(y, AAtrx);
  }
}

#endif

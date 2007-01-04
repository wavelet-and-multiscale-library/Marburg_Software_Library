// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MATRIX_BLOCK_H
#define _MATHTL_MATRIX_BLOCK_H

#include <iostream>

namespace MathTL
{
  /*!
    Abstract base class for all those matrices which shall be used as
    matrix blocks in the class BlockMatrix.
  */
  template <class C>
  class MatrixBlock
  {
  public:
    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef typename Vector<C>::size_type size_type;

    //! purely virtual destructor
    virtual ~MatrixBlock() = 0;

    //! row dimension
    virtual const size_type row_dimension() const = 0;
    
    //! column dimension
    virtual const size_type column_dimension() const = 0;

    /*!
      matrix-vector multiplication Mx = (*this) * x;
      we assume that the vector Mx has the correct size and
      is not identical to x
    */
    virtual void apply(const Vector<C>& x, Vector<C>& Mx) const = 0;
    
    /*!
      transposed matrix-vector multiplication Mtx = (*this)^T * x;
      we assume that the vector Mtx has the correct size and
      is not identical to x
    */
    virtual void apply_transposed(const Vector<C>& x, Vector<C>& Mtx) const = 0;

    //! clone the matrix block
    virtual MatrixBlock<C>* clone() const = 0;

    //! get a transposed version of the matrix block
    virtual MatrixBlock<C>* clone_transposed() const = 0;

    //! print block onto a stream
    virtual void print(std::ostream& os,
		       const unsigned int tabwidth = 10,
		       const unsigned int precision = 3) const = 0;
  };
}

#include <algebra/matrix_block.cpp>

#endif

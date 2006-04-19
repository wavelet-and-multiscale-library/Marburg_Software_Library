// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
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

    //! print block onto a stream
    virtual void print(std::ostream& os,
		       const unsigned int tabwidth = 10,
		       const unsigned int precision = 3) const = 0;
  };
}

#include <algebra/matrix_block.cpp>

#endif

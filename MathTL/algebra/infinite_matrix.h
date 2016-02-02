// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_INFINITE_MATRIX_H
#define _MATHTL_INFINITE_MATRIX_H

namespace MathTL
{
  /*!
    abstract base class for infinite-dimensional diagonal matrices
    (that can serve for scaling an InfiniteVector<C,I>...)
  */
  template <class C, class I = int>
  class InfiniteDiagonalMatrix
  {
  public:
    /*!
      destructor
    */
    virtual ~InfiniteDiagonalMatrix() = 0;

    /*!
      evaluate lambda-th entry on the diagonal
    */
    virtual double diag(const I& lambda) const = 0;
  };
}

#include <algebra/infinite_matrix.cpp>

#endif

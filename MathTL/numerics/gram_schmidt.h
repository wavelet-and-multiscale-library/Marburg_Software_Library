// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_GRAM_SCHMIDT_H
#define _MATHTL_GRAM_SCHMIDT_H

#include<utils/array1d.h>

namespace MathTL
{
  //! Gram-Schmidt orthogonalization process
  /*!
    Gram-Schmidt process, orthogonalize a set of linear independent
    vectors using projection

    The "vector" class V ist supposed to have the member functions
    add, scale and inner_product.
   */
  template <class V>
  void gramSchmidtProcess(Array1D<V>& vectors);
}

#include <numerics/gram_schmidt.cpp>

#endif // _MATHTL_GRAM_SCHMIDT_H

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_DIFFERENCES_H
#define _MATHTL_DIFFERENCES_H

#include <iostream>
#include <algebra/infinite_vector.h>

namespace MathTL
{
  /*!
    Apply a k-th forward difference to an infinite vector of function values
    (which we assume to correspond to a uniformly sampled mapping on the integers).
  */
  template <unsigned int K>
  InfiniteVector<double, int> forward_difference(const InfiniteVector<double, int>& a);

  /*!
    the 0-th forward difference is the identity
   */
  template <>
  InfiniteVector<double, int> forward_difference<0>(const InfiniteVector<double, int>& a);

  /*!
    first forward difference
  */
  template <>
  InfiniteVector<double, int> forward_difference<1>(const InfiniteVector<double, int>& a);

  /*!
    Apply a k-th backward difference to an infinite vector of function values.
  */
  template <unsigned int K>
  InfiniteVector<double, int> backward_difference(const InfiniteVector<double, int>& a);

  /*!
    the 0-th backward difference is the identity
  */
  template <>
  InfiniteVector<double, int> backward_difference<0>(const InfiniteVector<double, int>& a);

  /*!
    first backward difference
  */
  template <>
  InfiniteVector<double, int> backward_difference<1>(const InfiniteVector<double, int>& a);
}

// include implementation
#include <numerics/differences.cpp>

#endif

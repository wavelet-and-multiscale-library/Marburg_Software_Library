// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MULTI_DIFFERENCES_H
#define _MATHTL_MULTI_DIFFERENCES_H

#include <iostream>
#include <algebra/infinite_vector.h>
#include <utils/multiindex.h>

namespace MathTL
{
  /*!
    Apply a first forward difference in the nu-th direction
    to an infinite vector of function values
    (which we assume to correspond to a uniformly sampled mapping on the integers).
  */
  template <unsigned int DIMENSION, unsigned int DIRECTION>
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  multivariate_forward_difference(const InfiniteVector<double, MultiIndex<int, DIMENSION> >& a);

  /*!
    Apply a first backward difference in the nu-th direction
    to an infinite vector of function values.
  */
  template <unsigned int DIMENSION, unsigned int DIRECTION>
  InfiniteVector<double, MultiIndex<int, DIMENSION> >
  multivariate_backward_difference(const InfiniteVector<double, MultiIndex<int, DIMENSION> >& a);
}

// include implementation
#include <numerics/multi_differences.cpp>

#endif

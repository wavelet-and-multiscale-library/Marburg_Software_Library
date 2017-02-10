// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_STEEPEST_DESCENT1_H
#define _FRAME_TL_STEEPEST_DESCENT1_H

#include <algebra/infinite_vector.h>
#include <algebra/vector.h>

namespace FrameTL
{
  /*!
  */
  template <class PROBLEM>
  void steepest_descent1_SOLVE(const PROBLEM& P, const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& u_epsilon, 
			      const int jmax, InfiniteVector<double,typename PROBLEM::Index> rhs = InfiniteVector<double,typename PROBLEM::Index>());
}

#include <steepest_descent1.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_STEEPEST_DESCENT_H
#define _FRAME_TL_STEEPEST_DESCENT_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{
  /*!
  */
  template <class PROBLEM>
  void steepest_descent_SOLVE(const PROBLEM& P, const double epsilon,
			      InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);
}

#include <steepest_descent.cpp>

#endif

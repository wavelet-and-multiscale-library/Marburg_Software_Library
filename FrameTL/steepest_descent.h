// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2010                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_STEEPEST_DESCENT_H
#define _FRAME_TL_STEEPEST_DESCENT_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{
  
  /*!
    \file steepest_descent.h
    The adaptive steepest descent wavelet frame algorithm from Dahlke,
    Fornasier, Raasch, Stevenson, Werner 2007
   */

  /*!
    \brief
    Adaptive steepest descent wavelet frame algorithm from Dahlke, Fornasier, Raasch,
    Stevenson, Werner 2007.
    \param P The cached discrete problem.
    \param epsilon The target \f$\ell_2\f$-accuracy of the algorithm.
    \param approximations An array of length number of patches+1. We return in this array
    the local discrete approximations on each patch.
    The last entry contains the final global discrete approximation at termination.
  */
  template <class PROBLEM>
  void steepest_descent_SOLVE(const PROBLEM& P, const double epsilon,
			      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);
}

#include <steepest_descent.cpp>

#endif

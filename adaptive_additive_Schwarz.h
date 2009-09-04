// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2010                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_ADDITIVE_H
#define _FRAME_TL_ADDITIVE_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{

  /*!
    \file adaptive_additive_Schwarz.h
    The adaptive additive Schwarz solver.
  */
  /*!
    \brief
    Adaptive additive Schwarz wavelet frame algorithm from PhD thesis Werner 2009.
    \param P The cached discrete problem.
    \param epsilon The target \f$\ell_2\f$-accuracy of the algorithm.
    \param approximations An array of length number of patches+1. We return in this array
    the local discrete approximations on each patch.}
    The last entry contains the final global discrete approximation at termination.
   */
  template <class PROBLEM>
  void AddSchw(const PROBLEM& P, const double epsilon,
	       Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);


}

#include <adaptive_additive_Schwarz.cpp>

#endif

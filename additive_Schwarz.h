// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_ADDITIVE_H
#define _FRAME_TL_ADDITTVE_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{
  /*!
  */
  template <class PROBLEM>
  void addtive_Schwarz_SOLVE(const PROBLEM& P, const double epsilon,
				    InfiniteVector<double, typename PROBLEM::Index>& u_epsilon);
}

#include <additive_Schwarz.cpp>

#endif

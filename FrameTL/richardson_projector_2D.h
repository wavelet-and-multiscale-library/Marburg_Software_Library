// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_RICHARDSON_PROJECTOR_2D_H
#define _FRAME_TL_RICHARDSON_PROJECTOR_2D_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{


  /*!
  */
  template <class PROBLEM, class HILFSPROBLEM>
  void richardson_SOLVE_projector_2D(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);

 

}

#include <richardson_projector_2D.cpp>

#endif

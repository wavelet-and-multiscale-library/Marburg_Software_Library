// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_MULTIPLIKATIVE_H
#define _FRAME_TL_MULTIPLIKATIVE_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{


  /*!
    adaptive multiplicative Schwarz method from [DRSW]
   */
  template <class PROBLEM>
  void  MultSchw(const PROBLEM& P, const double epsilon,
		 Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);


  /*!
    adaptive multiplicative Schwarz method from [DRSW]
    tuned version
   */
  template <class PROBLEM>
  void  adaptive_multiplicative_Schwarz_SOLVE(const PROBLEM& P,  const double epsilon,
					      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);

}

#include <adaptive_multiplicative_Schwarz.cpp>

#endif

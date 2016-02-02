// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2010                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAME_TL_MULTIPLIKATIVE_H
#define _FRAME_TL_MULTIPLIKATIVE_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{
  /*!
    \file adaptive_multiplicative_Schwarz.h
    The adaptive multiplicative Schwarz solver.
  */
  
  /*!
    \brief  Adaptive multiplicative Schwarz wavelet frame algorithm from Stevenson, Werner 2009.

    \param P The cached discrete problem.
    \param epsilon The target \f$\ell_2\f$-accuracy of the algorithm.
    \param approximations An array of length number of patches +1. We return in this array
    the local discrete approximations on each patch.
    The last entry contains the final global discrete approximation at termination.
  */
  template <class PROBLEM>
  void  MultSchw_Proj(const PROBLEM& P, const double epsilon,
		 Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);

//   /*!
//     adaptive multiplicative Schwarz method from [DRSW]
//     tuned version
//    */
//   template <class PROBLEM>
//   void  adaptive_multiplicative_Schwarz_SOLVE(const PROBLEM& P,  const double epsilon,
// 					      Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);

}

#include <adaptive_multiplicative_Schwarz_Proj.cpp>

#endif

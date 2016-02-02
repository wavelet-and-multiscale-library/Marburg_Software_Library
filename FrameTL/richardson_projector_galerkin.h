// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+
//
// edited: Dominik Lellek, 2010
//
// an implementable version of MODSOLVE (Stevenson) 
// local systems are solved using the Galerkin method from CDD

#ifndef _FRAME_TL_RICHARDSON_PROJECTOR_GALERKIN_H
#define _FRAME_TL_RICHARDSON_PROJECTOR_GALERKIN_H

#include <algebra/infinite_vector.h>

namespace FrameTL
{


  /*!
  */
  template <class PROBLEM>
  void richardson_SOLVE_projector_galerkin(const PROBLEM& P, const double epsilon,
			InfiniteVector<double, typename PROBLEM::Index>& u_epsilon,
			Array1D<InfiniteVector<double, typename PROBLEM::Index> >& approximations);

 

}

#include <richardson_projector_galerkin.cpp>

#endif

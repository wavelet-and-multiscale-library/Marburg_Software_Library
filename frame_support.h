// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_FRAME_SUPPORT_H
#define _FRAMETL_FRAME_SUPPORT_H

#include <list>
#include <set>

//#include <geometry/chart.h>
#include <aggregated_frame.h>

namespace FrameTL
{

  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AffineLinearMapping<DIM_d>& ch1,
			  const AffineLinearMapping<DIM_d>&ch2,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu);
 
  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const LinearBezierMapping& ch1,
			  const LinearBezierMapping& ch2,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu);
}
#include <frame_support.cpp>
  
#endif

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

#include <geometry/chart.h>
#include <geometry/point.h>
#include <cube/cube_basis.h>

using WaveletTL::CubeBasis;

namespace FrameTL
{


  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const Point<DIM_m>& p);

  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  typename CubeBasis<IBASIS,DIM_d>::Support& supp_lambda,
		  const Point<DIM_m>& p);

  /*!
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  typename CubeBasis<IBASIS,DIM_d>::Support& supp_lambda,
			  typename CubeBasis<IBASIS,DIM_d>::Support& supp_mu);

  /*!
    tests wether the line segments defined by the points
    A and B as well as C and D intersect.
     0 = no intersection
     1 = infinitely many intersection points
     2 = single intersection point, but at least one knot involved
     3 = single intersection point situated 
         in the inner of both line segments
   */
  template <unsigned int DIM>
  int edgesIntersect (const Point<DIM>& A, const Point<DIM>& B,
		      const Point<DIM>& C, const Point<DIM>& D);


  /*!
    Checks whether p lies left of right of or on the line specified by the
    Points p1 and p2. This line's orientation is given by the vector starting in p1 and
    ending in p2.
    returning 0 means RIGHT OF LINE
    returning 1 means LEFT OF LINE
    returning 2 means ON LINE
   */
  template <unsigned int DIM>
  unsigned short int pos_wrt_line (const Point<DIM>& p,
				   const Point<DIM>& p1, const Point<DIM>&  p2);
  
 
}
#include <frame_support.cpp>
  
#endif

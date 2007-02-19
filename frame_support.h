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
#include <galerkin/infinite_preconditioner.h>
#include <index1D.h>



using WaveletTL::CubeBasis;
using FrameTL::Index1D;

namespace FrameTL
{
  /*!
    checks wether Point 'p' lies in the support of
    the wavelet \psi_\lambda
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const Point<DIM_m>& p);

  /*!
    checks wether Point 'p' lies in the support of
    the wavelet \psi_\lambda, when the support of the
    corresponding cube wavelet is already known and
    thus the latter is passed as an argument
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
		  const Point<DIM_m>& p);

  /*!
    checks wether the support of the wavelet frame elements intersect
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu);

  /*!
    checks wether the support of the wavelet frame elements intersect   
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda);


  /*!
    THIS ROUTINE IS INTENDED FOR THE SPECIAL CASE OF TRIVIAL
    PARAMETRIZATIONS, NAMELY AFFINE LINEAR MAPPINGS 'A x + B' WITH 'A' BEEING
    A DIAGONAL MATRIX.
    The function checks wether two wavelets intersect and returns
    an irregualar partition of the support intersection pulled back
    to the unit cube by the chart corresponding to \psi_\lambda.
    This is needed to be able to exactly compute the entries
    of the stiffness matrix for the above case of very simple patch
    parametrizations.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu,
			  FixedArray1D<Array1D<double>,DIM_d >& supp_intersect);

  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports_1D(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
 			     const Index1D<IBASIS>& lambda,
 			     const Index1D<IBASIS>& mu,
 			     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
 			     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu,
			     const int dir,
			     Array1D<double>& supp_intersect);
  

  /*!
    For a given wavelet frame element \psi_\lambda, compute all generators/wavelets
    \psi_\nu with level |\nu|=j, such that the respective supports
    have a nontrivial intersection
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void intersecting_wavelets(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting);

  /*!
    Decide whether the support of a given (primal) generator/wavelet \psi_\lambda
    intersects the singular support of another (primal) generator/wavelet \psi_\nu.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_singular_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& nu);



  template <unsigned int DIM>
  bool qudarangles_intersect (FixedArray1D<Point<DIM>, 4> poly1, FixedArray1D<Point<DIM>, 4> poly2);

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
  unsigned short int pos_wrt_line (const Point<DIM>& p, const Point<DIM>& p1, const Point<DIM>&  p2);
  
 
}
#include <frame_support.cpp>
  
#endif

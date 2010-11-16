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
#include <aggregated_frame.h>

using WaveletTL::CubeBasis;
using FrameTL::Index1D;

namespace FrameTL
{

//        Attention the calculation of the support between two Patches only works if the cube is mapped on a rectangle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  /*!
    \file frame_support.h
    Routines for handling the supports of aggregated wavelet frame elements.
  */


//   //
//   //  Checks wether Point 'p' lies in the support of the wavelet \f$\psi_\lambda\f$,
//   //  when the support of the corresponding cube wavelet is already known and
//   //  the latter is passed as an argument.
//   //
//   template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
//   bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
// 		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
// 		  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
// 		  const Point<DIM_m>& p);
  
   /*!
     Checks wether Point 'p' lies in the support of
     the wavelet frame element \f$\psi_\lambda\f$.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool in_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
		  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
		  const Point<DIM_m>& p);
  
  /*!
    Checks wether the support of the wavelet frame elements \f$\psi_\lambda\f$
    and \f$\psi_\mu\f$ intersect.
    The supports of the corresponding cube wavelets have to be passed as arguments
    supp_lambda and supp_mu.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu);

  /*!
    For the case of rectangular patches, the supports of all wavelets frame elements
    between minimal and maximal level are computed. This routine is supposed to be called during
    the initialization process of an AggregatedFrame in its constructor.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void precompute_supports_simple(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
				  Array1D<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Support>& all_patch_supports);
  /*!
    Checks wether the support of the wavelet frame elements intersect.
    WE ASSUME THAT THE PATCHES ARE RECTANGULAR AND THAT THEY ARE ALIGNED WITH THE
    CARTESIAN GRID.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports_simple(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				 const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				 const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu);

  /*!
    Gives the support of a wavelet or generator.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				 const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				 typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Support& SuppType);

  /*!
    Gives the support of a wavelet or generator on the unit cube.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void support_on_cube(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				 const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				 typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Support& SuppType);


  /*!
    Checks wether the support of the wavelet frame elements intersect.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda);

  /*!
    This function checks whether the supports of the two functions given by lambda and mu
    intersect. It also generates an irregular partition of the support intersection pulled back
    to the unit interval by the parametric mapping corresponding to lambda.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports_1D(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
 			     const Index1D<IBASIS>& lambda,
 			     const Index1D<IBASIS>& mu,
 			     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
 			     const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu,
			     const int dir,
			     Array1D<double>& supp_intersect);

  /*!
    THIS ROUTINE IS INTENDED FOR THE SPECIAL CASE OF TRIVIAL
    PARAMETRIZATIONS, NAMELY AFFINE LINEAR MAPPINGS 'A x + B' WITH 'A' BEEING
    A DIAGONAL MATRIX.
    The function checks whether two wavelets intersect and returns an irregular partition of the
    support intersection pulled back to the unit cube by the chart corresponding to \f$\psi_\lambda\f$.
    This is needed to be able to exactly compute the entries of the stiffness matrix for the
    above case of very simple patch parametrizations.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_supports(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& mu,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_lambda,
			  const typename CubeBasis<IBASIS,DIM_d>::Support* supp_mu,
			  FixedArray1D<Array1D<double>,DIM_d >& supp_intersect);

  /*!
    For a given wavelet frame element \f$\psi_\lambda\f$, compute all generators/wavelets
    \f$\psi_\nu\f$ with level \f$|\nu|=j\f$, such that the respective supports
    have a nontrivial intersection.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void intersecting_wavelets(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			     const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting);

  /*!
    For a given wavelet frame element \f$\psi_\lambda\f$, compute all generators/wavelets
    \f$\psi_\nu\f$ with level \f$|\nu|=j\f$ on patch p, such that the respective supports
    have a nontrivial intersection.
    WE ASSUME AGAIN THAT ONLY RECTANGULAR PATCHTES ALIGNED WITH THE COORDINATE AXES ARE USED.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  void intersecting_wavelets_on_patch(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				      const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				      const int p,
				      const int j, const bool generators,
				      std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting);

  /*!
    For a given wavelet frame element \f$\psi_\lambda\f$, compute all generators/wavelets from patch p
    in the set Lambda, the supports of which intersect the one of \f$\psi_\lambda\f$.
    WE ASSUME AGAIN THAT ONLY RECTANGULAR PATCHTES ALIGNED WITH THE COORDINATE AXES ARE USED.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  inline
  void intersecting_wavelets (const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			      const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			      const int p,
			      const std::set<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& Lambda,	      
			      std::list<typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& intersecting);

  /*!
    Decide whether the support of a given generator/wavelet \f$\psi_\lambda\f$
    intersects the singular support of another generator/wavelet \f$\psi_\nu\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  bool intersect_singular_support(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
				  const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& nu);

  /*!
    This routine tests whether the line segments defined by the points
    A and B as well as C and D intersect.
     0 = no intersection
     1 = infinitely many intersection points
     2 = single intersection point, but at least one knot involved
     3 = single intersection point situated in the inner of both line segments
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

  /*!
    This function checks whether the convex qudrangles given by the vertices in poly1 and poly2
    have a non-trivial intersection.
   */
  template <unsigned int DIM>
  bool quadrangles_intersect (FixedArray1D<Point<DIM>, 4> poly1, FixedArray1D<Point<DIM>, 4> poly2);
  
 
}
#include <frame_support.cpp>
  
#endif

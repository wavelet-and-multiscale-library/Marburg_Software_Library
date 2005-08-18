// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_IVP_H
#define _MATHTL_IVP_H

#include <geometry/point.h>
#include <algebra/matrix.h>

namespace MathTL
{
  /*!
    abstract base class for a general (vector valued)
    two-point boundary value problem on [0,1]
    
      y'(t) = f(t,y(t))

    with linear, homogeneous boundary conditions

      r(y(0),y(1)) = Ay(0) + By(1) = 0
  */
  template <unsigned int DIM>
  class BVP
  {
  public:
    /*!
      virtual destructor
     */
    virtual ~BVP ();

    /*!
      apply the left boundary condition matrix to v
    */
    virtual void apply_A(const Point<DIM>& v, Point<DIM>& result) const = 0;
    
    /*!
      apply the right boundary condition matrix to v
    */
    virtual void apply_B(const Point<DIM>& v, Point<DIM>& result) const = 0;
    
    /*!
      apply the right--hand side f to (t,v)
    */
    virtual void apply_f(const double t, const Point<DIM>& v, Point<DIM>& result) const = 0;
  };
}

#include <numerics/bvp.cpp>

#endif

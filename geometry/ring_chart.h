// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_RING_CHART_H
#define _MATHTL_RING_CHART_H

#include <geometry/chart.h>

using std::string;

namespace MathTL
{
  /*!
    A parametrization of the ring-shaped domain
      R = {(x,y) : r_0 <= ||(x,y)|| <= r_1 }
    by
      kappa: (0,1)^2 -> R,
      kappa(s,phi) = r(s)(cos(2*pi*phi),sin(2*pi*phi))
    with
      r(s) = r_0+s*(r_1-r_0).
   */
  class RingChart
    : public Chart<2,2>
  {
  public:
    //! constructor from two radii
    RingChart(const double r0, const double r1);

    //! virtual destructor
    virtual ~RingChart();

    //! map a point (forward) y = kappa(x)
    void map_point(const Point<2>& x, Point<2>& y) const;

    //! dummy 1D relict
    double map_point(const double, const int) const { return 0; }

    //! inverse mapping y = kappa^{-1}(x)
    void map_point_inv(const Point<2>& x, Point<2>& y) const;

    //! dummy 1D relict
    double map_point_inv(const double, const int) const { return 0; }
    
    /*!
      square root of the Gram determinant sqrt(det(Dkappa(x)^T * Dkappa(x)))
      (additional factor for integration over "plain" functions)
    */
    const double Gram_factor(const Point<2>& x) const;

    /*!
      i-th partial derivative of Gram factor
      (additional factor for integration over first derivatives)
    */
    const double Gram_D_factor(const unsigned int i,
			       const Point<2>& x) const;

    //! (i,j)-th element of (D (kappa^{-1}))(x)
    const double Dkappa_inv(const unsigned int i,
				    const unsigned int j,
				    const Point<2>& x) const;

    //! checks whether a special point x lies in the patch represented by this parametrization
    const bool in_patch(const Point<2>& x) const;

    //! dummy
    const double a_i(const int i) const { return 1.0; };    

    //! returns a string representation of this object
    const string to_string() const;
    
  protected:
    double r0_, r1_;
  };
  
}

#include "geometry/ring_chart.cpp"

#endif

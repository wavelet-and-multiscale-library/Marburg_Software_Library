// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_CHART_H
#define _MATHTL_CHART_H

#include <iostream>
#include <algebra/matrix.h>
#include <algebra/vector.h>
#include <geometry/point.h>

namespace MathTL
{
  /*!
    Abstract base class for smooth parametrizations
      kappa: (0,1)^d -> R^m
    of single "patches" in R^m.
   */
  template <unsigned int DIM_d, unsigned int DIM_m>
  class Chart
  {
  public:
    //! virtual destructor
    virtual ~Chart() {}

    /*!
      map a point (forward) y = kappa(x)
     */
    virtual void map_point(const Point<DIM_d>& x, Point<DIM_m>& y) const = 0;

    /*!
      inverse mapping y = kappa^{-1}(x)
     */
    virtual void map_point_inv(const Point<DIM_d>& x, Point<DIM_m>& y) const = 0;

    /*!
      square root of the Gram determinant sqrt(det(Dkappa(x)^T * Dkappa(x)))
      (additional factor for integration over "plain" functions)
    */
    virtual const double Gram_factor(const Point<DIM_d>& x) const = 0;

    /*!
      i-th partial derivative of Gram factor
      (additional factor for integration over first derivatives)
    */
    virtual const double Gram_D_factor(const unsigned int i,
				       const Point<DIM_d>& x) const = 0;

    /*!
      (i,j)-th element of the inverse of Dkappa at x (Dkappa)^{-1}(x)
    */
    virtual const double Dkappa_inv(const unsigned int i,
				    const unsigned int j,
				    const Point<DIM_d>& x) const = 0;

    /*!
      checks whether a special point x lies in the patch represented by this
      parametrization
    */
    virtual const bool in_patch(const Point<DIM_m>& x) const = 0;
  };
  
  //
  // Some examples:

  //! identity mapping
  template <unsigned int DIM>
  class IdentityMapping
    : public Chart<DIM,DIM>
  {
  public:
    void map_point(const Point<DIM>&, Point<DIM>&) const;
    void map_point_inv(const Point<DIM>&, Point<DIM>&) const;
    const double Gram_factor(const Point<DIM>&) const;
    const double Gram_D_factor(const unsigned int i, const Point<DIM>& x) const;
    const double Dkappa_inv(const unsigned int i, const unsigned int j,
			    const Point<DIM>& x) const;
    const bool in_patch(const Point<DIM>& x) const;
  };

  /*!
    affine linear mapping y = A*x+b
   */
  template <unsigned int DIM>
  class AffineLinearMapping
    : public Chart<DIM,DIM>
  {
  public:
    //! constructor from A and b (dimensions should fit)
    AffineLinearMapping(const Matrix<double>& A, const Point<DIM>& b);

    void map_point(const Point<DIM>&, Point<DIM>&) const;
    void map_point_inv(const Point<DIM>&, Point<DIM>&) const;
    const double Gram_factor(const Point<DIM>&) const;
    const double Gram_D_factor(const unsigned int i, const Point<DIM>& x) const;
    const double Dkappa_inv(const unsigned int i, const unsigned int j,
			    const Point<DIM>& x) const;
    const bool in_patch(const Point<DIM>& x) const;

    //! read access to A
    const Matrix<double>& A() const { return A_; }
    
    //! read access to b
    const Point<DIM>& b() const { return b_; }

  protected:
    Matrix<double> A_, A_inv;
    double det_A;
    Point<DIM> b_;
  };
}

#include "geometry/chart.cpp"

#endif

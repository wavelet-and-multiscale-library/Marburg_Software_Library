// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SPLINES_H
#define _MATHTL_SPLINES_H

#include <cmath>
#include <iostream>
#include <utils/array1d.h>
#include <utils/function.h>

namespace MathTL
{
  /*!
    This class models (compactly supported) splines of order d
      f(x) = sum_{j=0}^n alpha_j B_{j,d}(x)
    with respect to the nondecreasing knot sequence
      t_0 <= t_1 <= ... <= t_n,
    where
      B_{j,d}(x) = (t_{j+d}-t_j)[t_j,...,t_{j+d}](.-x)^{d-1}_+
    is the j-th (normalized) B-spline of order d.

    References:
    * deBoor: A practical guide to splines
  */
  template <unsigned int d>
  class Spline : public Function<1>
  {
  public:
    /*!
      default constructor: splines are real-valued;
      the default knot sequence is {0,1,...,d} to house a cardinal B-spline
    */
    Spline()
      : Function<1>(1), knots_(d+1)
    {
      for (unsigned int i(0); i <= d; i++)
	knots_[i] = i;
    }

    /*!
      construct spline from a knot sequence
     */
    Spline(const Array1D<double>& knots)
      : Function<1>(1), knots_(knots)
    {
    }

    /*!
      virtual destructor
    */
    virtual ~Spline() {}

    /*!
      value of a spline
    */
    double value(const Point<1>& p,
		 const unsigned int component = 0) const
    {
      return 0;
    }
  
    /*!
      value of a spline
    */
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }

    void dump()
    {
      std::cout << knots_ << std::endl;
    }
    
  protected:
    /*!
      the knot sequence t_0,...,t_n
     */
    Array1D<double> knots_;
  };
}

#endif

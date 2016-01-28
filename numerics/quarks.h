// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_QUARKS_H
#define _MATHTL_QUARKS_H

#include <cmath>
#include <utils/function.h>
#include <numerics/splines.h>
#include <utils/tiny_tools.h>
#include <algebra/piecewise.h>
#include <algebra/polynomial.h>
#include <algebra/sparse_matrix.h>
#include <numerics/cardinal_splines.h>

namespace MathTL
{
  /*!
    evaluate a shifted and dilatated Quark-function
      phi_{p,j,k}(x) = ((2^jx-k)/((d+1)/2))^p2^{j/2}N_d(2^jx-k+d/2)
  */
  template <int d>
  inline
  double EvaluateQuark_td(const int p, const int j, const int k, const double x){
  return pow(((pow(2,j)*x-k)/((d+1)/2)),p)*EvaluateCardinalBSpline_td<d>(j, k, x);
  }


template <int d, int p, int j, int k>
  class Quark_td : public Function<1>
  {
  public:
    /*!
      default constructor: B-spline-quarks are real-valued
    */
    Quark_td() : Function<1>(1) {}

    /*!
      virtual destructor
    */
    virtual ~Quark_td() {}

    /*!
      value of a B-spline-quark
    */
    inline double value(const Point<1>& x,
			const unsigned int component = 0) const
    {
      return EvaluateQuark_td<d>(p, j, k, x(0));
    }
  
    /*!
      value of a B-spline-quark
    */
    void vector_value(const Point<1> &x,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(x);
    }
  };
}
#endif

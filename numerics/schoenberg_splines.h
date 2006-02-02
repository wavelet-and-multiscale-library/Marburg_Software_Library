// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SCHOENBERG_SPLINES_H
#define _MATHTL_SCHOENBERG_SPLINES_H

#include <cmath>
#include <utils/function.h>

namespace MathTL
{
  /*
    The following routines provide the pointwise evaluation of
    d-th order Schoenberg B-splines N_{k,d}, k >= -d+1,
    which use the special infinite knot sequence

      t^0_{-d+1} = ... = t^0_0 = 0          (knot with multiplicity d at x=0)
      t^0_k = k, k >= 1

    i.e., t^0_k = max(0,k) for k >= -d+1.
  */

  /*!
    evaluate an arbitrary Schoenberg B-spline N_{k,d}(x)
  */
  template <int d>
  double EvaluateSchoenbergBSpline(const int k, const double x)
  {
    double r(0);
    
    // take care of the multiple knot t_{-d+1} = ... = t_0
    double diff = std::max(0,k+d-1) - std::max(0,k);
    if (diff > 0) r += (x - std::max(0,k)) * EvaluateSchoenbergBSpline<d-1>(k, x) / diff;
    diff = std::max(0,k+d) - std::max(0,k+1);
    if (diff > 0) r += (std::max(0,k+d) - x) * EvaluateSchoenbergBSpline<d-1>(k+1, x) / diff;
    
    return r;
  }
  
  /*!
    evaluate an arbitrary Schoenberg B-spline N_{k,1}(x) = \chi_{[t^0_k,t^0_{k+1})}
  */
  template <>
  inline
  double EvaluateSchoenbergBSpline<1>(const int k, const double x)
  {
    return (x >= std::max(0,k) && x < std::max(0,k+1) ? 1.0 : 0.0);
  }
  
  /*!
    evaluate the first derivative N_{k,d}'(x) of an arbitrary Schoenberg B-spline
  */
  template <int d>
  inline
  double EvaluateSchoenbergBSpline_x(const int k, const double x)
  {
    double r(0);

    if (k == 1-d) {
      r = (1-d) * EvaluateSchoenbergBSpline<d-1>(2-d, x);
    } else {
      if (k >= 0)
	r = EvaluateSchoenbergBSpline<d-1>(k, x) - EvaluateSchoenbergBSpline<d-1>(k+1, x);
      else
	r = (d-1) * (EvaluateSchoenbergBSpline<d-1>(k, x) / (k+d-1)
		     - EvaluateSchoenbergBSpline<d-1>(k+1, x) / (k+d));
    }

    return r;
  }
  
  /*!
    evaluate the first derivative N_{k,1}'(x) of an arbitrary Schoenberg B-spline
  */
  template <>
  inline
  double EvaluateSchoenbergBSpline_x<1>(const int k, const double x)
  {
    return 0.;
  }

}

#endif

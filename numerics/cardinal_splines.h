// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_CARDINAL_SPLINES_H
#define _MATHTL_CARDINAL_SPLINES_H

#include <cmath>
#include <utils/function.h>

namespace MathTL
{
  /*!
    evaluate a shifted cardinal B-spline N_d(x-k) via recursion
  */
  template <int d>
  double EvaluateCardinalBSpline(const int k, const double x)
  {
    return ((x-k) * EvaluateCardinalBSpline<d-1>(k, x)
	    + (k+d-x) * EvaluateCardinalBSpline<d-1>(k+1, x)) / (d-1);
  }

  /*!
    evaluate a shifted cardinal B-spline N_1(x-k)
  */
  template <>
  inline
  double EvaluateCardinalBSpline<1>(const int k, const double x)
  {
    if (x < k)
      return 0.;
    else
      if (x >= k+1)
	return 0.;
    return 1.;
  }

  /*!
    evaluate a primal CDF function
      phi_{j,k}(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_td(const int j, const int k, const double x)
  {
    const double factor(ldexp(1.0, j));
    return sqrt(factor) * EvaluateCardinalBSpline<d>(k, factor * x + d/2);
  }
  
  /*!
    evaluate the first derivative N_d'(x-k) of a shifted cardinal B-spline
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_x(const int k, const double x)
  {
    if (d == 1)
      return 0.;
    else
      return EvaluateCardinalBSpline<d-1>(k, x) - EvaluateCardinalBSpline<d-1>(k+1, x);
  }
  
  /*!
    evaluate the first derivative of a primal CDF function
      phi_{j,k}'(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  template <int d>
  inline
  double EvaluateCardinalBSpline_td_x(const int j, const int k, const double x)
  {
    const double factor(ldexp(1.0, j));
    return factor * sqrt(factor) * EvaluateCardinalBSpline_x<d>(k, factor * x + d/2);
  }

  /*!
    evaluate a shifted cardinal B-spline N_d(x-k) via recursion
    (remark: only use this version if the spline order d is unknown at compile time)
  */
  double EvaluateCardinalBSpline(const int d, const int k, const double x)
  {
    if (x < k)
      return 0.;
    else
      {
	if (x >= k+d)
	  return 0.;
	else
	  {
	    /* we know that x\in\supp N_d(.-k) */
	    if (d == 1)
	      return 1.;
	    else
	      {
		/* hard-encode case d=2 */
		if (d == 2)
		  {
		    if (x < k+1)
		      return x-k;
		    else
		      return 2.-(x-k);
		  }
		else
		  return ((x-k) * EvaluateCardinalBSpline(d-1, k, x)
			  + (k+d-x) * EvaluateCardinalBSpline(d-1, k+1, x)) / (d-1);
	      }
	  }
      }
  }
  
  /*!
    evaluate a primal CDF function
      phi_{j,k}(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  inline double EvaluateCardinalBSpline_td(const int d, const int j, const int k, const double x)
  {
    const double factor(ldexp(1.0, j));
    return sqrt(factor) * EvaluateCardinalBSpline(d, k, factor * x + d/2);
  }
  
  /*!
    evaluate the first derivative N_d'(x-k) of a shifted cardinal B-spline
  */
  inline double EvaluateCardinalBSpline_x(const int d, const int k, const double x)
  {
    if (d == 1)
      return 0.;
    else
      return EvaluateCardinalBSpline(d-1, k, x) - EvaluateCardinalBSpline(d-1, k+1, x);
  }
  
  /*!
    evaluate the first derivative of a primal CDF function
      phi_{j,k}'(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  inline double EvaluateCardinalBSpline_td_x(const int d, const int j, const int k, const double x)
  {
    const double factor(ldexp(1.0, j));
    return factor * sqrt(factor) * EvaluateCardinalBSpline_x(d, k, factor * x + d/2);
  }

  /*!
    cardinal B-spline N_d(x) as Function object
  */
  template <int d>
  class CardinalBSpline : public Function<1>
  {
  public:
    /*!
      default constructor: B-splines are real-valued
    */
    CardinalBSpline() : Function<1>(1) {}

    /*!
      virtual destructor
    */
    virtual ~CardinalBSpline() {}

    /*!
      value of a B-spline
    */
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      return EvaluateCardinalBSpline<d>(0, p(0));
    }
  
    /*!
      value of a B-spline
    */
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
  };
}

#endif

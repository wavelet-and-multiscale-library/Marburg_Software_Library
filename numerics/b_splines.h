// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_B_SPLINES_H
#define _MATHTL_B_SPLINES_H

#include <cmath>
#include <utils/function.h>

namespace MathTL
{
  /*!
    evaluate a cardinal B-spline N_d(x) via recursion
  */
  template <int d>
  double evaluate_Bspline(const double x)
  {
    return (x*evaluate_Bspline<d-1>(x)
      + (d-x)*evaluate_Bspline<d-1>(x-1))/(d-1);
  }
 
  /*!
    evaluate a cardinal B-spline N_1(x) = \chi_{[0,1)}
  */
  template <>
  inline
  double evaluate_Bspline<1>(const double x)
  {
    return (x >= 0 && x < 1.0 ? 1.0 : 0.0);
  }

  /*!
    evaluate a shifted cardinal B-spline N_d(x-k) via recursion
  */
  double EvaluateBspline(const int d, const int k, const double x)
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
		  return ((x-k) * EvaluateBspline(d-1, k, x)
			  + (k+d-x) * EvaluateBspline(d-1, k+1, x)) / (d-1);
	      }
	  }
      }
  }
  
  /*!
    evaluate a primal CDF function
      phi_{j,k}(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  inline double EvaluateBspline_td(const int d, const int j, const int k, const double x)
  {
    const double factor(ldexp(1.0, j));
    return sqrt(factor) * EvaluateBspline(d, k, factor * x + d/2);
  }
  
  /*!
    evaluate the first derivative N_d'(x-k) of a shifted cardinal B-spline
  */
  inline double EvaluateBspline_x(const int d, const int k, const double x)
  {
    if (d == 1)
      return 0.;
    else
      return EvaluateBspline(d-1, k, x) - EvaluateBspline(d-1, k+1, x);
  }
  
  /*!
    evaluate the first derivative of a primal CDF function
      phi_{j,k}'(x) = 2^{j/2}N_d(2^jx-k+d/2)
  */
  inline double EvaluateBspline_td_x(const int d, const int j, const int k, const double x)
  {
    const double factor(ldexp(1.0, j));
    return factor * sqrt(factor) * EvaluateBspline_x(d, k, factor * x + d/2);
  }

  /*!
    cardinal B-spline N_d(x) as Function object
  */
  template <int d>
  class Bspline : public Function<1>
  {
  public:
    /*!
      default constructor: B-splines are real-valued
    */
    Bspline() : Function<1>(1) {}

    /*!
      virtual destructor
    */
    virtual ~Bspline() {}

    /*!
      value of a B-spline
    */
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      return evaluate_Bspline<d>(p(0));
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

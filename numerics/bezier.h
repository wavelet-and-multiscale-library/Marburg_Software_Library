// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_BEZIER_H
#define _MATHTL_BEZIER_H

#include <cmath>

namespace MathTL
{
  /*!
    point evaluation of the k-th Bernstein basis polynomial of degree d

      B_{k,d}(x) = binomial{d}{k} * x^k * (1-x)^{d-k}

    at x, 0<=k<=d.
  */
  template <int d>
  double EvaluateBernsteinPolynomial(const int k, const double x)
  {
    return
      (1-x) * EvaluateBernsteinPolynomial<d-1>(k,  x)
      + x   * EvaluateBernsteinPolynomial<d-1>(k-1,x);
  }

  /*!
    evaluation of B_{k,0}(x)
  */
  template <>
  double EvaluateBernsteinPolynomial<0>(const int k, const double x)
  {
    return (k == 0 ? 1.0 : 0.0);
  }

  /*!
    point evaluation of the first derivative of B_{k,d}(x)
  */
  template <int d>
  double EvaluateBernsteinPolynomial_x(const int k, const double x)
  {
    return d * (EvaluateBernsteinPolynomial<d-1>(k-1, x)
		- EvaluateBernsteinPolynomial<d-1>(k, x));
  }

  /*!
    point evaluation of a quadratic polynomial p(x) in Bernstein representation
     
      p(x) = sum_{k=0}^2 B_{k,2}(x) * b_k
  */
  double EvaluateBernsteinPolynomial(const double b0, const double b1, const double b2,
				     const double x)
  {
    const double b10 = (1-x)*b0+x*b1;
    const double b11 = (1-x)*b1+x*b2;
    
    return (1-x)*b10+x*b11;
  }
  
  /*!
    point evaluation of a cubic polynomial p(x) in Bernstein representation
     
      p(x) = sum_{k=0}^3 B_{k,3}(x) * b_k
  */
  double EvaluateBernsteinPolynomial(const double b0, const double b1,
				     const double b2, const double b3,
				     const double x)
  {
    const double b10 = (1-x)*b0+x*b1;
    const double b11 = (1-x)*b1+x*b2;
    const double b12 = (1-x)*b2+x*b3;
    
    const double b20 = (1-x)*b10+x*b11;
    const double b21 = (1-x)*b11+x*b12;

    return (1-x)*b20+x*b21;
  }

  /*!
    point evaluation of the first derivative of a cubic polynomial p(x) in Bernstein representation
     
      p'(x) = sum_{k=0}^3 B'_{k,3}(x) * b_k
            = 3 * sum_{k=0}^3 (B_{k-1,2}(x)-B_{k,2}) * b_k
            = 3 * sum_{k=0}^2 B_{k,2}(x) * (b_{k+1}-b_k)
  */
  inline
  double EvaluateBernsteinPolynomial_x(const double b0, const double b1,
				       const double b2, const double b3,
				       const double x)
  {
    return EvaluateBernsteinPolynomial(b1-b0, b2-b1, b3-b2, x);
  }

  /*!
    point evaluation of the i-th cubic Hermite basic interpolatory spline
      i=0 : s(0)=1, s'(0)=0
      i!=0: s(0)=0, s'(0)=1
  */
  double EvaluateHermiteSpline(const int i, const double x) {
    if (i == 0) {
      // phi_0, Bezier coefficients are {0,0,1,1,1,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinPolynomial(0, 0, 1, 1, x+1.0)
		: EvaluateBernsteinPolynomial(1, 1, 0, 0, x));
      }
    } else {
      // phi_1, Bezier coefficients are {0,0,-1/3,0,1/3,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinPolynomial(0, 0, -1./3., 0, x+1.0)
		: EvaluateBernsteinPolynomial(0, 1./3., 0, 0, x));
      }
    }
  }

  /*!
    evaluate a translated and dilated version of the i-th cubic Hermite interpolant
  */
  double EvaluateHermiteSpline_td(const int i, const int j, const int k, const double x) {
    const double factor(ldexp(1.0, j));
    return sqrt(factor) * EvaluateHermiteSpline(i, factor * x - k);
  }
  
  /*!
    point evaluation of the first derivative of the i-th cubic Hermite basic interpolatory spline
      i=0 : s(0)=1, s'(0)=0
      i!=0: s(0)=0, s'(0)=1
  */
  double EvaluateHermiteSpline_x(const int i, const double x) {
    if (i == 0) {
      // phi_0, Bezier coefficients are {0,0,1,1,1,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinPolynomial_x(0, 0, 1, 1, x+1.0)
		: EvaluateBernsteinPolynomial_x(1, 1, 0, 0, x));
      }
    } else {
      // phi_1, Bezier coefficients are {0,0,-1/3,0,1/3,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinPolynomial_x(0, 0, -1./3., 0, x+1.0)
		: EvaluateBernsteinPolynomial_x(0, 1./3., 0, 0, x));
      }
    }
  }

  /*!
    evaluate the first derivative of a translated and dilated version of the i-th cubic Hermite interpolant
  */
  double EvaluateHermiteSpline_td_x(const int i, const int j, const int k, const double x) {
    const double factor(ldexp(1.0, j));
    return factor * sqrt(factor) * EvaluateHermiteSpline_x(i, factor * x - k);
  }
  
}

#include <numerics/bezier.cpp>

#endif

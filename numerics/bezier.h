// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_BEZIER_H
#define _MATHTL_BEZIER_H

#include <cmath>
#include <utils/function.h>
#include <utils/tiny_tools.h>

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
    return ((k < 0 || k > d)
	    ? 0
	    : (1-x) * EvaluateBernsteinPolynomial<d-1>(k,  x)
	    + x   * EvaluateBernsteinPolynomial<d-1>(k-1,x));
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
    point evaluation of the second derivative of B_{k,d}(x)
  */
  template <int d>
  double EvaluateBernsteinPolynomial_xx(const int k, const double x)
  {
    return d*(d-1) * (EvaluateBernsteinPolynomial<d-2>(k-2, x)
		      - 2*EvaluateBernsteinPolynomial<d-2>(k-1, x)
		      + EvaluateBernsteinPolynomial<d-2>(k, x));
  }

  /*!
    point evaluation of a quadratic polynomial p(x) in Bernstein representation
     
      p(x) = sum_{k=0}^2 B_{k,2}(x) * b_k
  */
  inline
  double EvaluateBernsteinRepresentation(const double b0, const double b1, const double b2,
					 const double x)
  {
    const double b10 = (1.-x)*b0+x*b1;
    const double b11 = (1.-x)*b1+x*b2;
    
    return (1.-x)*b10+x*b11;
  }
  
  /*!
    point evaluation of a cubic polynomial p(x) in Bernstein representation
     
      p(x) = sum_{k=0}^3 B_{k,3}(x) * b_k
  */
  inline
  double EvaluateBernsteinRepresentation(const double b0, const double b1,
					 const double b2, const double b3,
					 const double x)
  {
    const double b10 (b0+x*(b1-b0)); // = (1.-x)*b0+x*b1;
    const double b11 (b1+x*(b2-b1)); // = (1.-x)*b1+x*b2;
    const double b12 (b2+x*(b3-b2)); // = (1.-x)*b2+x*b3;
    
    const double b20 (b10+x*(b11-b10)); // = (1.-x)*b10+x*b11;
    const double b21 (b11+x*(b12-b11)); // = (1.-x)*b11+x*b12;

    return b20+x*(b21-b20); // = (1.-x)*b20+x*b21;
  }

  /*!
    point evaluation of the first derivative of a cubic polynomial p(x) in Bernstein representation
     
      p'(x) = sum_{k=0}^3 B'_{k,3}(x) * b_k
            = 3 * sum_{k=0}^3 (B_{k-1,2}(x)-B_{k,2}(x)) * b_k
            = 3 * sum_{k=0}^2 B_{k,2}(x) * (b_{k+1}-b_k)
  */
  inline
  double EvaluateBernsteinRepresentation_x(const double b0, const double b1,
					   const double b2, const double b3,
					   const double x)
  {
    return 3*EvaluateBernsteinRepresentation(b1-b0, b2-b1, b3-b2, x);
  }

  /*!
    point evaluation of the second derivative of a cubic polynomial p(x) in Bernstein representation
     
      p''(x) = sum_{k=0}^3 B''_{k,3}(x) * b_k
             = 6 * sum_{k=0}^3 (B_{k-2,1}(x)-2*B_{k_1,1}(x)+B_{k,1}(x)) * b_k
             = 6 * sum_{k=0}^1 B_{k,1}(x) * (b_{k+2}-2*b_{k+1}-b_k)
  */
  inline
  double EvaluateBernsteinRepresentation_xx(const double b0, const double b1,
					    const double b2, const double b3,
					    const double x)
  {
    return 6*((1-x)*(b2-2*b1+b0) + x*(b3-2*b2+b1));
  }
  
  /*!
    point evaluation of the i-th cubic Hermite basic interpolatory spline
      i=0 : s(0)=1, s'(0)=0
      i!=0: s(0)=0, s'(0)=1
  */
  double EvaluateHermiteSpline(const int i, const double x) {
    if (i == 0) {
      // phi_0, Bezier coefficients are {0,0,1,1,1,0,0}
      if (x <= -1. || x >= 1.) {
	return 0.;
      } else {
	return (x <= 0.
		? EvaluateBernsteinRepresentation(0., 0., 1., 1., x+1.)
// 		(x+1)*(x+1)*(1-2*x)
		: EvaluateBernsteinRepresentation(1., 1., 0., 0., x)
// 		(1-x)*(1-x)*(2*x+1)
		);
      }
    } else {
      // phi_1, Bezier coefficients are {0,0,-1/3,0,1/3,0,0}
      if (x <= -1. || x >= 1.) {
	return 0.;
      } else {
	return (x <= 0.
		? EvaluateBernsteinRepresentation(0., 0., -1./3., 0., x+1.)
		: EvaluateBernsteinRepresentation(0., 1./3., 0., 0., x));
      }
    }
  }

  /*!
    evaluate a translated and dilated version of the i-th cubic Hermite interpolant
  */
  inline
  double EvaluateHermiteSpline_td(const int i, const int j, const int k, const double x) {
#if 0
    const double factor(1<<j);
    return sqrt(factor) * EvaluateHermiteSpline(i, factor * x - k);
#else
    return twotothejhalf(j) * EvaluateHermiteSpline(i, (1<<j) * x - k);
#endif
  }
  
  /*!
    point evaluation of the first derivative of the i-th cubic Hermite basic interpolatory spline
      i=0 : s(0)=1, s'(0)=0
      i!=0: s(0)=0, s'(0)=1
  */
  inline
  double EvaluateHermiteSpline_x(const int i, const double x) {
    if (i == 0) {
      // phi_0, Bezier coefficients are {0,0,1,1,1,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinRepresentation_x(0, 0, 1, 1, x+1.0)
		: EvaluateBernsteinRepresentation_x(1, 1, 0, 0, x));
      }
    } else {
      // phi_1, Bezier coefficients are {0,0,-1/3,0,1/3,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinRepresentation_x(0, 0, -1./3., 0, x+1.0)
		: EvaluateBernsteinRepresentation_x(0, 1./3., 0, 0, x));
      }
    }
  }

  /*!
    point evaluation of the second derivative of the i-th cubic Hermite basic interpolatory spline
      i=0 : s(0)=1, s'(0)=0
      i!=0: s(0)=0, s'(0)=1
  */
  inline
  double EvaluateHermiteSpline_xx(const int i, const double x) {
    if (i == 0) {
      // phi_0, Bezier coefficients are {0,0,1,1,1,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinRepresentation_xx(0, 0, 1, 1, x+1.0)
		: EvaluateBernsteinRepresentation_xx(1, 1, 0, 0, x));
      }
    } else {
      // phi_1, Bezier coefficients are {0,0,-1/3,0,1/3,0,0}
      if (x <= -1 || x >= 1) {
	return 0;
      } else {
	return (x <= 0
		? EvaluateBernsteinRepresentation_xx(0, 0, -1./3., 0, x+1.0)
		: EvaluateBernsteinRepresentation_xx(0, 1./3., 0, 0, x));
      }
    }
  }

  /*!
    evaluate the first derivative of a translated and dilated version of the i-th cubic Hermite interpolant
  */
  inline
  double EvaluateHermiteSpline_td_x(const int i, const int j, const int k, const double x) {
#if 0
    return twotothejhalf(3*j) * EvaluateHermiteSpline_x(i, (1<<j) * x - k);
#else
    return (double)(1<<j) * (double)twotothejhalf(j) * EvaluateHermiteSpline_x(i, (1<<j) * x - k);
#endif
  }
  
  /*!
    evaluate the second derivative of a translated and dilated version of the i-th cubic Hermite interpolant
  */
  inline
  double EvaluateHermiteSpline_td_xx(const int i, const int j, const int k, const double x) {
#if 0
    return twotothejhalf(5*j) * EvaluateHermiteSpline_xx(i, (1<<j) * x - k);
#else
    const double factor = (double)(1<<j);
    return factor * factor * twotothejhalf(j) * EvaluateHermiteSpline_xx(i, (1<<j) * x - k);
#endif
  }

  /*!
    translated and dilated Hermite interpolant as a function object
   */
  class CubicHermiteInterpolant_td : public Function<1>
  {
  public:
    //! constructor from j, c, k
    CubicHermiteInterpolant_td(const int j, const int c, const int k)
      : j_(j), c_(c), k_(k) {}
    
    //! virtual destructor
    virtual ~CubicHermiteInterpolant_td() {}
    
    //! point value
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      return EvaluateHermiteSpline_td(c_, j_, k_, p[0]);
    }
    
    //! point value
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
    
    //! set j
    void set_j(const int j) { j_ = j; }
    
    //! set c
    void set_c(const int c) { c_ = c; }

    //! set k
    void set_k(const int k) { k_ = k; }
    
  protected:
    int j_, c_, k_;
  };
  
  /*!
    tensor product of translated and dilated Hermite interpolants as a function object
   */
  class CubicHermiteInterpolant2D_td : public Function<2>
  {
  public:
    //! constructor from j, c, k
    CubicHermiteInterpolant2D_td(const int j, const int c0, const int c1, const int k0, const int k1)
      : j_(j), c0_(c0), c1_(c1), k0_(k0), k1_(k1) {}
    
    //! virtual destructor
    virtual ~CubicHermiteInterpolant2D_td() {}
    
    //! point value
    inline double value(const Point<2>& p,
			const unsigned int component = 0) const
    {
      return EvaluateHermiteSpline_td(c0_, j_, k0_, p[0])
	* EvaluateHermiteSpline_td(c1_, j_, k1_, p[1]);
    }
    
    //! point value
    void vector_value(const Point<2> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
    
    //! set j
    void set_j(const int j) { j_ = j; }
    
    //! set c0
    void set_c0(const int c0) { c0_ = c0; }

    //! set c1
    void set_c1(const int c1) { c1_ = c1; }

    //! set k0
    void set_k0(const int k0) { k0_ = k0; }
    
    //! set k1
    void set_k1(const int k1) { k1_ = k1; }
    
  protected:
    int j_, c0_, c1_, k0_, k1_;
  };
  
}

#include <numerics/bezier.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_SPLINES_H
#define _MATHTL_SPLINES_H

#include <cmath>
#include <iostream>
#include <utils/array1d.h>
#include <utils/function.h>
#include <algebra/matrix.h>

namespace MathTL
{
  /*!
    abstract base class for (semi)infinite knot sequences in Z
      t_{k0} <= ... <= t_k <= t_{k+1} <= ...
    for splines the support of which overlaps [0,\infty)
  */
  class KnotSequence
  {
  public:
    virtual ~KnotSequence() = 0;
    
    /*
      index of the first knot
     */
    virtual int k0() const = 0;
    
    /*!
      compute the k-th knot
    */
    virtual double knot(const int k) const = 0;
  };

  KnotSequence::~KnotSequence() {}

  /*!
    evaluate an arbitrary B-spline N_{j,d}(x) via recursion
  */
  template <int d>
  double evaluate_Bspline(const KnotSequence* knots, const int j, const double x)
  {
    double r(0);

    // take care of multiple knots
    double diff = knots->knot(j+d-1) - knots->knot(j);
    if (diff > 0) r += (x - knots->knot(j)) * evaluate_Bspline<d-1>(knots, j, x) / diff;
    diff = knots->knot(j+d) - knots->knot(j+1);
    if (diff > 0) r += (knots->knot(j+d) - x) * evaluate_Bspline<d-1>(knots, j+1, x) / diff;
    
    return r;
  }
  
  /*!
    evaluate an arbitrary B-spline N_{j,1}(x) = \chi_{[t_j,t_{j+1})}
  */
  template <>
  inline
  double evaluate_Bspline<1>(const KnotSequence* knots, const int j, const double x)
  {
    return (x >= knots->knot(j) && x < knots->knot(j+1) ? 1.0 : 0.0);
  }

  /*!
    Given a (semi)infinite knot sequence in Z
      t_{k0} <= ... <= t_k <= t_{k+1} <= ...
    with t_k = k for k >= 0, such that the knot sequence
    (t_j)_j is a subset of ((1/2)*t_j)_j,
    we know that the corresponding d-th order B-splines N_{k,d} fulfill a
    refinement relation of the form
      N_{k0+k,d}(x) = \sum_{n=0}^{d-2-k0} m_{n,k} N_{k0+n,d}(2x), 0<=k<=-1-k0

    This routine computes the refinement matrix M=(m_{n,k})_{n,k} for n=0..d-2-k0, k=0..-1-k0.

    References:
    [P] Primbs:
        Stabile biorthogonale Wavelet-Basen auf dem Intervall,
	Dissertation, Univ. Duisburg-Essen, 2006
  */
  template <int d>
  void
  compute_Bspline_refinement_matrix(const KnotSequence* knots, Matrix<double>& M);


  /*!
    evaluate an arbitrary B-spline N_{j,d}(x) via recursion
  */
  template <int d>
  double evaluate_Bspline(const Array1D<double>& knots, const unsigned int j, const double x)
  {
    double r(0);

    // take care of multiple knots
    double diff = knots[j+d-1] - knots[j];
    if (diff > 0) r += (x - knots[j]) * evaluate_Bspline<d-1>(knots, j, x) / diff;
    diff = knots[j+d] - knots[j+1];
    if (diff > 0) r += (knots[j+d] - x) * evaluate_Bspline<d-1>(knots, j+1, x) / diff;
    
    return r;
  }

  /*!
    evaluate an arbitrary B-spline N_{j,1}(x) = \chi_{[t_j,t_{j+1})}
  */
  template <>
  inline
  double evaluate_Bspline<1>(const Array1D<double>& knots, const unsigned int j, const double x)
  {
    return (x >= knots[j] && x < knots[j+1] ? 1.0 : 0.0);
  }

  /*!
    This class models (compactly supported) splines of order d
      f(x) = sum_{j=0}^n alpha_j N_{j,d}(x)
    with respect to the nondecreasing knot sequence
      t_0 <= t_1 <= ... <= t_{n+d},
    where
      N_{j,d}(x) = (t_{j+d}-t_j)[t_j,...,t_{j+d}](.-x)^{d-1}_+
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
      the default knot sequence is {0,1,...,d},
      the default spline is a cardinal B-spline, alpha_j = delta_{j,0}
    */
    Spline()
      : Function<1>(1), knots_(d+1), coeffs_(1)
    {
      for (unsigned int i(0); i <= d; i++)
	knots_[i] = i;

      coeffs_[0] = 1.0;
    }

    /*!
      construct spline from a knot sequence and given coefficients
     */
    Spline(const Array1D<double>& knots, const Array1D<double>& coeffs)
      : Function<1>(1), knots_(knots), coeffs_(coeffs)
    {
      assert(knots.size() == coeffs.size()+d);
    }

    /*!
      virtual destructor
    */
    virtual ~Spline() {}

    /*!
      value of a spline
    */
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      double r(0);
      
      for (unsigned int i(0); i < coeffs_.size(); i++) {
	if (coeffs_[i] != 0)
	  r += coeffs_[i] * evaluate_Bspline<d>(knots_, i, p(0));
      }

      return r;
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
    
  protected:
    /*!
      knot sequence t_0,...,t_{n+d}
     */
    Array1D<double> knots_;

    /*!
      B-spline coefficients
    */
    Array1D<double> coeffs_;
  };
}

#include <numerics/splines.cpp>

#endif

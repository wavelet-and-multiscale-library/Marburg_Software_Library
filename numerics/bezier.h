// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_BEZIER_H
#define _MATHTL_BEZIER_H

namespace MathTL
{
  /*!
    point evaluation of the k-th Bernstein polynomial of degree d

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
}

#include <numerics/bezier.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_GAUSS_QUADRATURE_H
#define _MATHTL_GAUSS_QUADRATURE_H

#include <numerics/quadrature.h>

namespace MathTL
{
  /*!
    N-point 1D Gauss-Legendre quadrature rule
    (uses precomputed points and weights for 1<=N<=10
  */
  template <unsigned int N = 1>
  class GaussLegendreRule
    : public QuadratureRule<1>
  {
  public:
    /*!
      construct N-point Gauss-legendre rule
    */
    GaussLegendreRule();
  };
}

// include implementation of inline functions
#include <numerics/gauss_quadrature.cpp>

#endif

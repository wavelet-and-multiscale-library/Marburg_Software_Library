// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MULTI_LAURENT_POLYNOMIAL_H
#define _MATHTL_MULTI_LAURENT_POLYNOMIAL_H

#include <cassert>
#include <utils/function.h>
#include <algebra/infinite_vector.h>
#include <algebra/polynomial.h>
#include <geometry/point.h>
#include <utils/multiindex.h>

namespace MathTL
{
  /*!
    A template class for general multivariate Laurent polynomials, i.e.,
    expressions of the form
      a(z)=\sum_{k\in\mathbb Z^d} a_k z^k
    The a_k are elements of a (number) ring R, z is from R^d.
  
    Examples: R=\mathbb R or \mathbb C
  */
  template <class R, unsigned int DIMENSION>
    class MultiLaurentPolynomial
    : protected InfiniteVector<R,MultiIndex<int, DIMENSION> >,
    public Function<DIMENSION, R>
    {
      
      /*!
	evaluate the Laurent polynomial (Horner scheme)
      */
      R value(const Point<DIMENSION>& p,
	      const unsigned int component = 0) const;
      
      /*!
	evaluate the Laurent polynomial (Horner scheme)
      */
      void vector_value(const Point<DIMENSION> &p,
			Vector<R>& values) const;
      
    };
}

// include implementation
#include <algebra/multi_laurent_polynomial.cpp>

#endif

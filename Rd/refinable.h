// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_REFINABLE_H
#define _WAVELETTL_REFINABLE_H

#include <iostream>
#include <geometry/sampled_mapping.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    Base class for a refinable function \phi on \mathbb R.

    A refinable function \phi on \mathbb R^d fulfills the refinement equation
    
      \phi(x) = \sum_{k\in\mathbb Z^d} a_k * \phi(M * x - k)

    where M\in\mathbb Z^{d\times d} is an expanding matrix and the a_k are
    real coefficients.
    In this module we will (at least for the moment) only support the case
    d = 1 and M = 2, so that

      \phi(x) = \sum_{k\in\mathbb Z} a_k * \phi(2 * x - k)

    The coefficients a_k are stored in a univariate Laurent polynomial.
    Since a refinable function admits much more functionality, we do not only
    use LaurentPolynomial<double> but specify the mask via a template parameter
    MASK which is assumed to have the signature of a LaurentPolynomial<double>.
  */
  template <class MASK>
  class RefinableFunction
    : public MASK
  {
  public:
    /*!
      Evaluate the n-th derivative of a dilated and translated version
      of the refinable function \phi
      
        (d/dx)^n \phi_{j,k}(x) = 2^{jn} * 2^{j/2} * \phi^{(n)}(2^j * x - k)
      
      on a dyadic subgrid of the interval [a,b].

      You should know at compile time which derivative is needed,
      since it is specified as a template paramter.

      We assume that \phi^{(n)} is zero at the boundary of its support.
    */
    template <unsigned int DERIVATIVE>
    SampledMapping<1> evaluate(const int j, const int k,
			       const int a, const int b,
			       const int resolution) const;

    /*!
      Compute the n-th moment
        M_n:=\int_{-\infty}^\infty x^n\phi(x)\,dx
      of the refinable function \phi
    */
    double moment(const unsigned int n) const;

  protected:
    /*!
      helper function for the moment calculation
    */
    double cnk(const unsigned int n, const unsigned int k) const;
  };
}

#include <Rd/refinable.cpp>

#endif

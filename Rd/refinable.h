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
#include <utils/multiindex.h>
#include <algebra/infinite_vector.h>
#include <geometry/point.h>
#include <geometry/sampled_mapping.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    Base class for a multivariate refinable function \phi on \mathbb R^d.

    A refinable function \phi on \mathbb R^d fulfills the refinement equation
    
      \phi(x) = \sum_{k\in\mathbb Z^d} a_k * \phi(M * x - k)

    where M\in\mathbb Z^{d\times d} is an expanding matrix and the a_k are
    real coefficients.
    In this module we will (at least for the moment) only support the case
    M = 2*I, so that

      \phi(x) = \sum_{k\in\mathbb Z^d} a_k * \phi(2 * x - k)

    The coefficients a_k are stored in a multivariate Laurent polynomial.
    Since a refinable function admits much more functionality, we do not only
    use MultivariateLaurentPolynomial<double,d> but specify the
    mask via a template parameter MASK which is assumed to have the signature
    of a MultivariateLaurentPolynomial<double, d>.
  */
  template <class MASK, unsigned int DIMENSION>
  class MultivariateRefinableFunction
    : public MASK
  {
  public:
    /*!
      Evaluate the refinable function \phi on the grid 2^{-resolution}\mathbb Z^d.
      We assume that \phi is zero at the boundary of its support.
    */
    InfiniteVector<double, MultiIndex<int, DIMENSION> >
    evaluate(const int resolution = 0) const;
    
    /*!
      Evaluate the mu-th (partial) derivative of the refinable function \phi
      on the grid 2^{-resolution}\mathbb Z^d.
      We assume that the derivative of \phi is zero at the boundary of its support (!)
    */
    InfiniteVector<double, MultiIndex<int, DIMENSION> >
    evaluate(const MultiIndex<unsigned int, DIMENSION>& mu,
	     const int resolution = 0) const;

    /*!
      Evaluate the mu-th (partial) derivative of a dilated and translated version
      of the refinable function \phi
      
      \partial^\mu \phi_{j,k}(x) = 2^{j|\mu|} * 2^{jd/2} * \partial^\mu\phi(2^j*x-k)
      
      on a dyadic subgrid of the cuboid [a,b]^d.

      We assume that \partial^\mu\phi is zero at the boundary of its support.
    */
    SampledMapping<DIMENSION>
    evaluate(const MultiIndex<unsigned int, DIMENSION>& mu,
	     const int j,
	     const MultiIndex<int, DIMENSION>& k,
	     const MultiIndex<int, DIMENSION>& a,
	     const MultiIndex<int, DIMENSION>& b,
	     const int resolution) const;

    /*!
      Compute the \alpha-th (monomial) moment
        M_\alpha:=\int_{-\infty}^\infty x^\alpha\phi(x)\,dx
      of the refinable function \phi

      Reference:
      [BBDK] Barinka, Barsch, Dahlke, Konik:
             Some Remarks on Quadrature Formulae for Refinable Functions and Wavelets
    */
    double moment(const MultiIndex<unsigned int, DIMENSION>& alpha) const;
  };

  /*!
    shorthand notation for univariate refinable functions
  */
  template <class MASK>
  class RefinableFunction
    : public MultivariateRefinableFunction<MASK, 1>
  {
  };

}

#include <Rd/refinable.cpp>

#endif

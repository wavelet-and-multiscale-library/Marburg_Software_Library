// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_GAUSS_QUADRATURE_H
#define _MATHTL_GAUSS_QUADRATURE_H

#include <utils/array1d.h>
#include <numerics/quadrature.h>
#include <numerics/ortho_poly.h>

namespace MathTL
{
  /*!
    N-point 1D Gauss-Legendre quadrature rule on [0,1]
    (uses precomputed points and weights for 1<=N<=10
  */
  class GaussLegendreRule
    : public QuadratureRule<1>
  {
  public:
    /*!
      construct N-point Gauss-legendre rule on [0,1]
    */
    GaussLegendreRule(const unsigned int N);
  };

  /*!
    N-point 1D Gauss rule for arbitrary weight functions,
    several constructors are available:
    - from the coefficients a_k, b_k in the three-term recursion
        P_k(x) = (x-alpha_k) * P_{k-1}(x) - beta_k * P_{k-2}(x),   k=1,2,...
      for the w-orthogonal polynomials P_k with leading coefficient 1,
      where p_{-1}(t)=0, p_0(t)=1
    - from 2N (generalized) moments of the form
        \nu_k = \int_a^b T_k(x)w(x)dx
      where the polynomials T_k also fulfill a recurrence relation
        T_k(x) = (x-a_k) * T_{k-1}(x) - b_k * T_{k-2}(x),   k=1,2,...

    references:
    * Sack/Donovan: An Algorithm for Gaussian Quadrature given Modified
      Moments, Numer. Math. 18(1972), 465-478
    * Golub/Gutknecht: Modified Moments for Indefinite Weight Functions,
      Numer. Math. 57(1990), 607-624
  */
  class GaussRule
    : public QuadratureRule<1>
  {
  public:
    /*!
      construct N-point Gauss rule from three-term recursion coefficients
      for orthogonal polynomials P on [a,b]

      \param P orth. polynomial, should provide a_1,...,a_N and b_1,...,b_N
      \param N number of quadrature points
    */
    GaussRule(const OrthogonalPolynomial& P,
	      const double a, const double b,
	      const unsigned int N);

    /*!
      construct N-point Gauss rule from (at least) 2N (monomial) moments
      of the weight function \int_a^b x^k w(x)dx
    */
    GaussRule(const Array1D<double>& moments,
	      const double a, const double b,
	      const unsigned int N);

    /*!
      construct N-point Gauss rule from (at least) 2N generalized moments
      of the weight function \int_a^b T_k(x)w(x)dx
    */
    GaussRule(const Array1D<double>& moments,
	      const OrthogonalPolynomial& T,
	      const double a, const double b,
	      const unsigned int N);

  private:
    /*!
      shared initialization routine
    */
    void init(const OrthogonalPolynomial& P,
	      const double a, const double b,
	      const unsigned int N);
  };
}

// include implementation of inline functions
#include <numerics/gauss_quadrature.cpp>

#endif

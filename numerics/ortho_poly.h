// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ORTHO_POLY_H
#define _MATHTL_ORTHO_POLY_H

#include <cmath>
#include <algebra/vector.h>
#include <algebra/polynomial.h>

namespace MathTL
{
  // some tools for polynomials that fulfill a three-term recurrence
  // relation (such as orthogonal polynomials, but also the monomials)
  //
  // reference: Deuflhard/Bornemann, Numerische Mathematik 1

  /*!
    Abstract base class for orthogonal polynomials that fulfill a
    (homogeneous) three-term recursion of the form

      p_k(t) = (t+a_k) * p_{k-1}(t) + b_k * p_{k-2}(t),   k=1,2,...

    where p_{-1}(t)=0, p_0(t)=1.
    Due to the specific shape of the recurrence relation, the
    p_k will have leading coefficient 1.
   */
  class OrthogonalPolynomial
  {
  public:
    /*!
      the coefficients a_k
    */
    virtual double a(const unsigned int k) const = 0;

    /*!
      the coefficients b_k
    */
    virtual double b(const unsigned int k) const = 0;

    /*!
      evaluate n-th orthogonal polynomial at x
    */
    double operator () (const unsigned int n, const double x) const;

    /*!
      assemble n-th orthogonal polynomial
    */
    Polynomial<double> assemble(const unsigned int n) const;

    /*!
      (trivial) forward summation of
        \sum_{k=0}^n \alpha_k * p_k(x)
     */
    double forwardSummation(const Vector<double>& coeffs, const double x) const;

    /*!
      adjoint summation of 
        \sum_{k=0}^n \alpha_k * p_k(x)
      remarks:
      - generalized Horner scheme (-> less multiplications than forward summation)
      - numerically stable for dominant solutions of the three-term recursion
      - potentially unstable for minimal solutions
    */
    double adjointSummation(const Vector<double>& coeffs, const double x) const;
  };

  //! monomials x^n
  class Monomial : public OrthogonalPolynomial
  {
  public:
    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
  };

  //! (normalized) Chebyshev polynomials
  class ChebyshevPolynomial : public OrthogonalPolynomial
  {
  public:
    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
  };

  //! (normalized) Legendre polynomials
  class LegendrePolynomial : public OrthogonalPolynomial
  {
  public:
    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
  };
}

// include implementation of inline functions
#include <numerics/ortho_poly.cpp>

#endif

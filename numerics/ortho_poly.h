// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ORTHO_POLY_H
#define _MATHTL_ORTHO_POLY_H

#include <cmath>
#include <utils/array1d.h>
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

      p_k(t) = (t-a_k) * p_{k-1}(t) - b_k * p_{k-2}(t),   k=1,2,...

    where p_{-1}(t)=0, p_0(t)=1.

    Due to the specific shape of the recurrence relation, the
    p_k will have leading coefficient 1. Moreover, we have

      a_k = <tp_{k-1}, p_{k-1}> / <p_{k-1}, p_{k-1}>
      b_k = <tp_{k-1}, p_{k-2}> / <p_{k-2}, p_{k-2}>
   */
  class OrthogonalPolynomial
  {
  public:
    //! virtual destructor
    virtual ~OrthogonalPolynomial() {}

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

  /*!
    the monomials x^n
  */
  class Monomial : public OrthogonalPolynomial
  {
  public:
    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
  };

  /*!
    Chebyshev polynomials on [-1,1] with leading coefficient 1
    (corresponding weight function: w(x)=1/sqrt(1-x^2) )
  */
  class ChebyshevPolynomial : public OrthogonalPolynomial
  {
  public:
    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
  };

  /*!
    Legendre polynomials on [-1,1] with leading coefficient 1
    (corresponding weight function: w(x)=1)
  */
  class LegendrePolynomial : public OrthogonalPolynomial
  {
  public:
    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
  };

  /*!
    orthogonal polynomials on [a,b] given by monomial or
    generalized moments of the weight function
      \nu_k = \int_a^b T_k(x)w(x)dx
    where the T_k also fulfill a three-term recursion;
    since we can only input a finite number of moments (2N),
    it is only possible to access a_1,...,a_N and b_1,...,b_N,
    so you should initialize with N as large as required for
    your purposes.

    references:
    + Sack/Donovan: An Algorithm for Gaussian Quadrature given Modified
      Moments, Numer. Math. 18(1972), 465-478
    + Golub/Gutknecht: Modified Moments for Indefinite Weight Functions,
      Numer. Math. 57(1990), 607-624
    + Golub/Welsch: Calculation of Gauss quadrature rules,
      Math. Comp. 23(1969), 221-230
  */
  class GenMomentsPolynomial : public OrthogonalPolynomial
  {
  public:
    /*!
      Golub/Welsch algorithm: [GW]
      
      \param moments (at least) 2*N+1 monomial moments of the weight function
      \param N number of three-term recursion coefficients created
    */
    GenMomentsPolynomial(const Array1D<double>& moments,
			 const double a, const double b,
			 const unsigned int N);
    
    /*!
      Sack/Donovan algorithm: 2(i),(ii) of [GG]
    */
    GenMomentsPolynomial(const Array1D<double>& moments,
			 const OrthogonalPolynomial& T,
			 const double a, const double b,
			 const unsigned int N);

    virtual ~GenMomentsPolynomial() {}

    double a(const unsigned int k) const;
    double b(const unsigned int k) const;
    
  private:
    /*!
      storage for the precomputed a_k and b_k
    */
    Array1D<double> ak_, bk_;
  };
}

// include implementation of inline functions
#include <numerics/ortho_poly.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_POLYNOMIAL_H
#define _MATHTL_POLYNOMIAL_H

#include <iostream>
#include <algebra/vector.h>
#include <utils/function.h>

namespace MathTL
{
  /*!
    A template class for univariate polynomials of the form
      p(x)=\sum_{k=0}^n a_k*x^k
    with coefficients from a class C.

    You can perform
    - point evaluation (via Horner scheme)
    - summation (coeffientwise)
    - multiplication with constants
    - multiplication with other polynomials
    - division with remainder
    - raising a polynomial to some power
    - substitution into another polynomial, special version for scale and shift
    - symbolic differentiation and integration
    - exact integration over an interval (with or without quadrature formulae)

    We derive Polynomial<C> from the class Function<1>, as it is indeed
    one. Since polynomials are completely determined by their coefficients,
    we derive the class also from Vector<C>, but protected to hide the Vector
    signature.
  */
  template <class C>
  class Polynomial
    : protected Vector<C>, public Function<1, C>
  {
  public:
    /*!
      default constructor: yields the zero polynomial of degree zero
    */
    Polynomial();

    /*!
      copy constructor
    */
    Polynomial(const Polynomial<C>& p);

    /*!
      constructor from a coefficient vector
    */
    Polynomial(const Vector<C>& coeffs);

    /*!
      constructor from a constant
    */
    explicit Polynomial(const C c);

    /*!
      virtual destructor
    */
    virtual ~Polynomial();

    /*!
      degree of the polynomial
    */
    unsigned int degree() const;

    /*!
      read-only access to single coefficients
    */
    C get_coefficient(const unsigned int k) const;

    /*!
      get all coefficients at once
    */
    void get_coefficients(Vector<C>& coeffs) const;

    /*!
      write access to single coefficients
    */
    void set_coefficient(const unsigned int k, const C coeff);

    /*!
      set all coefficients at once
    */
    void set_coefficients(const Vector<C>& coeffs);

    /*!
      evaluate the polynomial (Horner scheme)
    */
    C value(const C x) const;

    /*!
      evaluate n-th derivative of the polynomial (full Horner scheme)
    */
    C value(const C x, const unsigned int derivative) const;

    /*!
      evaluate the polynomial (Horner scheme)
      (calls the above value(const C))
    */
    C value(const Point<1>& p,
	    const unsigned int component = 0) const;

    /*!
      evaluate the polynomial (Horner scheme)
      (calls the above value(const C))
    */
    void vector_value(const Point<1> &p,
                      Vector<C>& values) const;

    /*!
      in place scaling p(x) -> p(s*x)
    */
    void scale(const C s);

    /*!
      in place shift p(x) -> p(x+s)
    */
    void shift(const C s);

    /*!
      substitute another polynomial into this one *this <- *this\circ p
    */
    void chain(const Polynomial<C>& p);

    /*!
      substitute another polynomial into this one without changing this object
      ( return *this\circ p )
    */
    Polynomial<C> substitute_into(const Polynomial<C>& p) const;

    /*!
      assignment of another polynomial
    */
    Polynomial<C>& operator = (const Polynomial<C>& p);

    /*!
      assignment of a constant polynomial (implicit conversion)
    */
    Polynomial<C>& operator = (const C c);

    /*!
      pointwise sum of two polynomials *this += p
    */
    void add(const Polynomial<C>& p);

    /*!
      pointwise weighted sum of two polynomials *this += s*p
    */
    void add(const C s, const Polynomial<C>& p);

    /*!
      pointwise weighted sum of two polynomials *this = s*(*this) + v
    */
    void sadd(const C s, const Polynomial<C>& p);

    /*!
      pointwise sum of two polynomials
    */
    Polynomial<C>& operator += (const Polynomial<C>& p);

    /*!
      pointwise sum of two polynomials
      (don't use this extensively, since one copy has to be made!)
    */
    Polynomial<C> operator + (const Polynomial<C>& p) const;

    /*!
      pointwise difference of two polynomials *this -= p
    */
    void subtract(const Polynomial<C> &p);

    /*!
      pointwise difference of two polynomials
    */
    Polynomial<C>& operator -= (const Polynomial<C> &p);

    /*!
      sign
      (makes a copy of *this)
    */
    Polynomial<C> operator - () const;
    
    /*!
      pointwise difference of two polynomials
      (don't use this extensively, since one copy has to be made!)
    */
    Polynomial<C> operator - (const Polynomial<C> &p) const;

    /*!
      multiplication with a real number
    */
    Polynomial<C>& operator *= (const C c);

    /*!
      multiplication with a real number
      (don't use this extensively, since one copy has to be made!)
    */
    Polynomial<C> operator * (const C c) const;

    /*!
      pointwise multiplication with another polynomial *this *= p
    */
    void multiply(const Polynomial<C>& p);

    /*!
      pointwise multiplication with another polynomial
    */
    Polynomial<C>& operator *= (const Polynomial<C>& p);

    /*!
      pointwise multiplication with another polynomial
      (don't use this extensively, since one copy has to be made!)
    */
    Polynomial<C> operator * (const Polynomial<C>& p);

    /*!
      raise the polynomial to some power
    */
    Polynomial<C> power(const unsigned int k) const;

    /*!
      divide the polynomial by another one with remainder: *this = p * q + r
    */
    void divide(const Polynomial<C>& q, Polynomial<C>& p) const;

    /*!
      divide the polynomial by another one with remainder: *this = p * q + r
    */
    void divide(const Polynomial<C>& q, Polynomial<C>& p, Polynomial<C>& r) const;

    /*!
      (symbolic) differentiation
    */
    Polynomial differentiate() const;

    /*!
      (symbolic) integration
    */
    Polynomial integrate() const;

    /*!
      integration over [a,b], optional use of (Gauss) quadrature formulae
    */
    double integrate(const double a,
                     const double b,
                     const bool quadrature = false) const;

    /*!
      inner product
    */
    double inner_product(Polynomial<C> p2, double a, double b);

  protected:
    /*!
      trim leading zero coefficients in the polynomial
    */
    void trim();
  };
  
  /*!
    multiplication with a real number
    (makes a copy)
  */
  template <class C>
  Polynomial<C> operator * (const C c, const Polynomial<C>& p);
  
  /*!
    stream output for polynomials,
    for readability we also print out the powers
  */
  template <class C>
  std::ostream& operator << (std::ostream &s, const Polynomial<C> &p);
}

// implementation of inline functions
#include <algebra/polynomial.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
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
    - raising a polynomial to some power
    - substitution into another polynomial
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
      write access to single coefficients
    */
    void set_coefficient(const unsigned int k, const C coeff);

    /*!
      evaluate the polynomial (Horner scheme)
    */
    C value(const Point<1>& p,
	    const unsigned int component = 0) const;

    /*!
      evaluate the polynomial (Horner scheme)
    */
    void vector_value(const Point<1> &p,
		      Array1D<C>& values) const;
  };

//   class Polynomial
//   {
//   public:
//     Polynomial(const Polynomial &p);
//     Polynomial(const CoeffsType *coeffs);
//     Polynomial(const double c);
//     ~Polynomial();
  
//     const CoeffsType *getCoeffs() const;
//     void              setCoeffs(const CoeffsType *coeffs);


//     Polynomial   sdifferentiate() const;                // symbolic differentiation
//     Polynomial   sintegrate() const;                    // symbolic integration
  
//     // integration over [a,b], optional use of quadrature formulae
//     double       integrate(const double a,
// 			   const double b,
// 			   const bool quadrature = false) const;
  
//     void         substituteIntoMe(const Polynomial &p);     // result: this\circ p
//     Polynomial   substituteInto(const Polynomial &p) const; // "

//     Polynomial power(const int k) const; // raise Polynomial to the k-th power

//     double operator () (const double x) const;  // point evaluation (Horner scheme)
//     double derivative(const double x) const;    // point evaluation of first derivative (Horner scheme)

//     Polynomial& operator = (const Polynomial &p);
//     Polynomial& operator = (const double c);                       // type conversion
  
//     Polynomial& operator += (const Polynomial &p); // pointwise
//     Polynomial& operator -= (const Polynomial &p); // pointwise
//     Polynomial& operator *= (const double c);      // multiplication with a constant
//     Polynomial& operator *= (const Polynomial &p); // pointwise

//     friend Polynomial operator + (const Polynomial &p, const Polynomial &q);
//     friend Polynomial operator - (const Polynomial &p, const Polynomial &q);
//     friend Polynomial operator - (const Polynomial &p);
//     friend Polynomial operator * (const double c, const Polynomial &p);  
//     friend Polynomial operator * (const Polynomial &p, const Polynomial &q);  
//   };

//   Polynomial operator + (const Polynomial &p, const Polynomial &q);
//   Polynomial operator - (const Polynomial &p, const Polynomial &q);
//   Polynomial operator - (const Polynomial &p);
//   Polynomial operator * (const double c, const Polynomial &p);  
//   Polynomial operator * (const Polynomial &p, const Polynomial &q);  

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

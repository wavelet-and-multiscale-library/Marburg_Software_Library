// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_LAURENT_POLYNOMIAL_H
#define _MATHTL_LAURENT_POLYNOMIAL_H

#include <cassert>
#include <utils/function.h>
#include <algebra/infinite_vector.h>
#include <algebra/polynomial.h>
#include <geometry/point.h>

namespace MathTL
{
  /*!
    A template class for general univariate Laurent polynomials, i.e.,
    expressions of the form
      a(z)=\sum_{k\in\mathbb Z} a_k z^k
    The a_k and z are elements of a (number) ring R.
  
    Examples: R=\mathbb R or \mathbb C
  */
  template <class R>
  class LaurentPolynomial
    : protected InfiniteVector<R,int>, public Function<1, R>
  {
  public:
    /*!
      const_iterator scanning the nontrivial coefficients
    */
    typedef typename InfiniteVector<R,int>::const_iterator
    const_iterator;

    /*!
      const_reverse_iterator scanning the nontrivial coefficients
    */
    typedef typename InfiniteVector<R,int>::const_reverse_iterator
    const_reverse_iterator;

    /*!
      default constructor, yields zero (Laurent) polynomial
    */
    LaurentPolynomial();

    /*!
      copy constructor
    */
    LaurentPolynomial(const LaurentPolynomial<R>& p);

    /*!
      constructor from a constant
    */
    explicit LaurentPolynomial(const R c);

    /*!
      constructor from a polynomial
    */
    explicit LaurentPolynomial(const Polynomial<R>& p);

    /*!
      virtual destructor
    */
    virtual ~LaurentPolynomial();

    /*!
      assignment of another Laurent polynomial
    */
    LaurentPolynomial<R>& operator = (const LaurentPolynomial<R>& p);

    /*!
      assignment of a constant
    */
    LaurentPolynomial<R>& operator = (const R c);

    /*!
      (Euclidean) degree of a Laurent polynomial
    */
    unsigned int degree() const;

    /*!
      (Polynomial-like) read-only access to single coefficients
    */
    R get_coefficient(const int k) const;

    /*!
      (Polynomial-like) write access to single coefficients
    */
    void set_coefficient(const int k, const R coeff);

    /*!
      const_iterator pointing to the first nontrivial coefficient
    */
    const_iterator begin() const;

    /*!
      const_iterator pointing to one behind the last nontrivial coefficient
    */
    const_iterator end() const;

    /*!
      const_reverse_iterator pointing to the last nontrivial coefficient
    */
    const_reverse_iterator rbegin() const;

    /*!
      const_reverse_iterator pointing to one before the first nontrivial coefficient
    */
    const_reverse_iterator rend() const;

    /*!
      evaluate the Laurent polynomial (Horner scheme)
    */
    R value(const R x) const;

    /*!
      evaluate the Laurent polynomial (Horner scheme)
      (calls the above value(const R))
    */
    R value(const Point<1>& p,
	    const unsigned int component = 0) const;

    /*!
      evaluate the Laurent polynomial (Horner scheme)
      (calls the above value(const R))
    */
    void vector_value(const Point<1> &p,
		      Vector<R>& values) const;

    /*!
      pointwise sum of two Laurent polynomials *this += p
    */
    void add(const LaurentPolynomial<R>& p);

    /*!
      pointwise weighted sum of two Laurent polynomials *this += s*p
    */
    void add(const R s, const LaurentPolynomial<R>& p);

    /*!
      pointwise weighted sum of two Laurent polynomials *this = s*(*this) + p
    */
    void sadd(const R s, const LaurentPolynomial<R>& p);

    /*!
      pointwise sum of two Laurent polynomials
    */
    LaurentPolynomial<R>& operator += (const LaurentPolynomial<R>& p);

    /*!
      pointwise sum of two Laurent polynomials
      (don't use this extensively, since one copy has to be made!)
    */
    LaurentPolynomial<R> operator + (const LaurentPolynomial<R>& p) const;

    /*!
      pointwise difference of two Laurent polynomials *this -= p
    */
    void subtract(const LaurentPolynomial<R>& p);

    /*!
      pointwise difference of two Laurent polynomials
    */
    LaurentPolynomial<R>& operator -= (const LaurentPolynomial<R> &p);

    /*!
      sign
      (makes a copy of *this)
    */
    LaurentPolynomial<R> operator - () const;
    
    /*!
      pointwise difference of two Laurent polynomials
      (don't use this extensively, since one copy has to be made!)
    */
    LaurentPolynomial<R> operator - (const LaurentPolynomial<R> &p) const;

    /*!
      multiplication with a real number
    */
    LaurentPolynomial<R>& operator *= (const R c);

    /*!
      multiplication with a real number from the right
      (don't use this extensively, since one copy has to be made!)
    */
    LaurentPolynomial<R> operator * (const R c) const;
    
    /*!
      pointwise multiplication with another Laurent polynomial
    */
    void multiply(const LaurentPolynomial<R>& p);

    /*!
      pointwise multiplication with another Laurent polynomial
    */
    LaurentPolynomial<R>& operator *= (const LaurentPolynomial<R>& p);

    /*!
      raise the Laurent polynomial to some power
    */
    LaurentPolynomial<R> power(const unsigned int k) const;

    /*!
      division with remainder by another Laurent polynomial q: *this = p * q + r
    */
    void divide(const LaurentPolynomial<R>& q, LaurentPolynomial<R>& p) const;

    /*!
      division with remainder by another Laurent polynomial q: *this = p * q + r
    */
    void divide(const LaurentPolynomial<R>& q,
		LaurentPolynomial<R>& p, LaurentPolynomial<R>& r) const;
  };

  /*!
    multiplication of a Laurent polynomial with a real number from the left
  */
  template <class R>
  LaurentPolynomial<R> operator * (const R c, const LaurentPolynomial<R>& p);

  /*!
    pointwise multiplication of two Laurent polynomials
  */
  template <class R>
  LaurentPolynomial<R> operator * (const LaurentPolynomial<R>& p,
				   const LaurentPolynomial<R>& q);
    
  /*!
    stream output for Laurent polynomials,
    for readability we also print out the powers
  */
  template <class R>
  std::ostream& operator << (std::ostream &s, const LaurentPolynomial<R> &p);

}

// include implementation of inline functions
#include <algebra/laurent_polynomial.cpp>

#endif

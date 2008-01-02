// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MULTI_LAURENT_POLYNOMIAL_H
#define _MATHTL_MULTI_LAURENT_POLYNOMIAL_H

#include <cassert>
#include <utils/function.h>
#include <algebra/infinite_vector.h>
#include <algebra/polynomial.h>
#include <algebra/laurent_polynomial.h>
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
  class MultivariateLaurentPolynomial
    : protected InfiniteVector<R, MultiIndex<int, DIMENSION> >,
      public Function<DIMENSION, R>
  {
  public:
    /*!
      const_iterator scanning the nontrivial coefficients
    */
    typedef typename InfiniteVector<R, MultiIndex<int, DIMENSION> >::const_iterator
    const_iterator;

    /*!
      const_reverse_iterator scanning the nontrivial coefficients
    */
    typedef typename InfiniteVector<R, MultiIndex<int, DIMENSION> >::const_reverse_iterator
    const_reverse_iterator;
    
    /*!
      default constructor, yields zero (Laurent) polynomial
    */
    MultivariateLaurentPolynomial();

    /*!
      copy constructor
    */
    MultivariateLaurentPolynomial(const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      constructor from a constant
    */
    explicit MultivariateLaurentPolynomial(const R c);

    /*!
      constructor from a univariate Laurent polynomial,
      embedding via z -> z_{n+1}
    */
    MultivariateLaurentPolynomial(const LaurentPolynomial<R>& p, const unsigned int n = 0);

    /*!
      virtual destructor
    */
    virtual ~MultivariateLaurentPolynomial();

    /*!
      assignment of another Laurent polynomial
    */
    MultivariateLaurentPolynomial<R, DIMENSION>& operator =
    (const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      assignment of a constant
    */
    MultivariateLaurentPolynomial<R, DIMENSION>& operator = (const R c);

    /*!
      (Polynomial-like) read-only access to single coefficients
    */
    R get_coefficient(const MultiIndex<int, DIMENSION>& k) const;

    /*!
      (Polynomial-like) write access to single coefficients
    */
    void set_coefficient(const MultiIndex<int, DIMENSION>& k, const R coeff);

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
      evaluate the Laurent polynomial
    */
    R value(const Point<DIMENSION>& p,
	    const unsigned int component = 0) const;
    
    /*!
      evaluate the Laurent polynomial
    */
    void vector_value(const Point<DIMENSION> &p,
		      Vector<R>& values) const;

    /*!
      pointwise sum of two Laurent polynomials *this += p
    */
    void add(const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      pointwise weighted sum of two Laurent polynomials *this += s*p
    */
    void add(const R s, const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      pointwise weighted sum of two Laurent polynomials *this = s*(*this) + p
    */
    void sadd(const R s, const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      pointwise sum of two Laurent polynomials
    */
    MultivariateLaurentPolynomial<R, DIMENSION>&
    operator += (const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      pointwise sum of two Laurent polynomials
      (don't use this extensively, since one copy has to be made!)
    */
    MultivariateLaurentPolynomial<R, DIMENSION>
    operator + (const MultivariateLaurentPolynomial<R, DIMENSION>& p) const;

    /*!
      pointwise difference of two Laurent polynomials *this -= p
    */
    void subtract(const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      pointwise difference of two Laurent polynomials
    */
    MultivariateLaurentPolynomial<R, DIMENSION>&
    operator -= (const MultivariateLaurentPolynomial<R, DIMENSION> &p);

    /*!
      sign
      (makes a copy of *this)
    */
    MultivariateLaurentPolynomial<R, DIMENSION> operator - () const;
    
    /*!
      pointwise difference of two Laurent polynomials
      (don't use this extensively, since one copy has to be made!)
    */
    MultivariateLaurentPolynomial<R, DIMENSION>
    operator - (const MultivariateLaurentPolynomial<R, DIMENSION> &p) const;

    /*!
      multiplication with a real number
    */
    MultivariateLaurentPolynomial<R, DIMENSION>& operator *= (const R c);

    /*!
      pointwise multiplication with another Laurent polynomial
    */
    void multiply(const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      pointwise multiplication with another Laurent polynomial
    */
    MultivariateLaurentPolynomial<R, DIMENSION>& operator *=
    (const MultivariateLaurentPolynomial<R, DIMENSION>& p);

    /*!
      just for testing
    */
    void dump() const;
  };

  /*!
    multiplication with a real number from the right
    (don't use this extensively, since one copy has to be made!)
  */
  template <class R, unsigned int DIMENSION>
  MultivariateLaurentPolynomial<R, DIMENSION> operator *
  (const R c, const MultivariateLaurentPolynomial<R, DIMENSION>& p);
  
  /*!
    pointwise multiplication of two Laurent polynomials
  */
  template <class R, unsigned int DIMENSION>
  MultivariateLaurentPolynomial<R, DIMENSION> operator *
  (const MultivariateLaurentPolynomial<R, DIMENSION>& p,
   const MultivariateLaurentPolynomial<R, DIMENSION>& q);

  /*!
    stream output for multivariate Laurent polynomials,
    for readability we also print out the powers
  */
  template <class R, unsigned int DIMENSION>
  std::ostream& operator <<
    (std::ostream& s, const MultivariateLaurentPolynomial<R, DIMENSION> &p);
}

// include implementation
#include <algebra/multi_laurent_polynomial.cpp>

#endif

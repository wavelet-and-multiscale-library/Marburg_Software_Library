// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_LAURENT_POLYNOMIAL_H
#define _MATHTL_LAURENT_POLYNOMIAL_H

#include <cassert>
#include <iterator>
#include <utils/function.h>
#include <algebra/infinite_vector.h>
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
    typedef typename InfiniteVector<R,int>::const_iterator const_iterator;

//     /*!
//       const_reverse_iterator scanning the nontrivial coefficients
//     */
//     typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;

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
      virtual destructor
    */
    virtual ~LaurentPolynomial();

    /*!
      assignment of another Laurent polynomial
    */
    LaurentPolynomial<R>& operator = (const LaurentPolynomial<R>& p);

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

//     /*!
//       const_reverse_iterator pointing to the last nontrivial coefficient
//     */
//     const_reverse_iterator rbegin() const;

//     /*!
//       const_reverse_iterator pointing to one before the first nontrivial coefficient
//     */
//     const_reverse_iterator rend() const;

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
      pointwise sum of two Laurent polynomials
    */
    LaurentPolynomial<R>& operator += (const LaurentPolynomial<R>& p);

    /*!
      pointwise sum of two Laurent polynomials
      (don't use this extensively, since one copy has to be made!)
    */
    LaurentPolynomial<R> operator + (const LaurentPolynomial<R>& p) const;

    /*!
      pointwise difference of two Laurentpolynomials
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
    
    
//     typedef typename std::map<int, R>::const_iterator LaurentIterator;
//     inline LaurentIterator begin() const { return _coeffs.begin(); }
//     inline LaurentIterator end() const { return _coeffs.end(); }

//     inline int kbegin() const { return _coeffs.begin()->first; }
//     inline int kend() const { return _coeffs.rbegin()->first; }

//     LaurentPolynomial<R>& operator -= (const LaurentPolynomial<R>& p);
//     LaurentPolynomial<R>& operator *= (const R c);
//     LaurentPolynomial<R>& operator *= (const LaurentPolynomial<R>& p);

//     friend LaurentPolynomial<R> operator + <> (const LaurentPolynomial<R>& p, const LaurentPolynomial<R>& q);
//     friend LaurentPolynomial<R> operator - <> (const LaurentPolynomial<R>& p, const LaurentPolynomial<R>& q);
//     friend LaurentPolynomial<R> operator - <> (const LaurentPolynomial<R>& p);
//     friend LaurentPolynomial<R> operator * <> (const R c, const LaurentPolynomial<R>& p);
//     friend LaurentPolynomial<R> operator * <> (const LaurentPolynomial<R>& p, const LaurentPolynomial<R>& q);

//     friend std::ostream& operator << <> (std::ostream& s, const LaurentPolynomial<R>&);
  
  private:
//     std::map<int, R> _coeffs;
  };

//   template <class R> LaurentPolynomial<R> power(const LaurentPolynomial<R>& p, const int n)
//   {
//     assert(n >= 0);

//     LaurentPolynomial<R> r(1);
//     for (int i(0); i < n; i++) r *= p;

//     return r;
//   }

  //
  //
  // implementation:


//   template <class R> LaurentPolynomial<R>& LaurentPolynomial<R>::operator = (const R& value)
//   {
//     _coeffs.clear();
//     setCoefficient(0, value);
//     return *this;
//   }


//   template <class R> LaurentPolynomial<R>& LaurentPolynomial<R>::operator += (const LaurentPolynomial<R>& p)
//   {
//     for (LaurentIterator it(p.begin()); it != p.end(); it++)
//       {
// 	R help(getCoefficient(it->first) + it->second);
// 	if (help != 0)
// 	  setCoefficient(it->first, help);
// 	else
// 	  _coeffs.erase(it->first);
//       }
//     return *this;
//   }

//   template <class R> LaurentPolynomial<R>& LaurentPolynomial<R>::operator -= (const LaurentPolynomial<R>& p)
//   {
//     for (LaurentIterator it(p.begin()); it != p.end(); it++)
//       {
// 	R help(getCoefficient(it->first) - it->second);
// 	if (help != 0)
// 	  setCoefficient(it->first, help);
// 	else
// 	  _coeffs.erase(it->first);
//       }
//     return *this;
//   }

//   template <class R> LaurentPolynomial<R>& LaurentPolynomial<R>::operator *= (const R c)
//   {
//     if (c == 0)
//       _coeffs.clear();
//     else
//       {
// 	typedef typename std::map<int, R>::iterator myIterator;
// 	for (myIterator it(_coeffs.begin()); it != _coeffs.end(); it++)
// 	  it->second *= c;
//       }
//     return *this;
//   }

//   template <class R> LaurentPolynomial<R>& LaurentPolynomial<R>::operator *= (const LaurentPolynomial<R>& p)
//   {
//     std::map<int, R> bak(_coeffs);

//     _coeffs.clear();
//     for (LaurentIterator itme(bak.begin()); itme != bak.end(); itme++)
//       for (LaurentIterator itp(p.begin()); itp != p.end(); itp++)
// 	{
// 	  _coeffs[itme->first+itp->first] += itme->second * itp->second;
// 	}

//     return *this;
//   }


//   //
//   //
//   // implementation of friend operators

//   template <class R> LaurentPolynomial<R> operator + (const LaurentPolynomial<R>& p, const LaurentPolynomial<R>& q)
//   {
//     LaurentPolynomial<R> r(p);
//     r += q;

//     return r;
//   }

//   template <class R> LaurentPolynomial<R> operator - (const LaurentPolynomial<R>& p, const LaurentPolynomial<R>& q)
//   {
//     LaurentPolynomial<R> r(p);
//     r -= q;

//     return r;
//   }

//   template <class R> LaurentPolynomial<R> operator - (const LaurentPolynomial<R>& p)
//   {
//     LaurentPolynomial<R> r(p);
//     r *= -1;

//     return r;
//   }

  /*!
    multiplication of a Laurent polynomial with a real number from the left
  */
  template <class R>
  inline
  LaurentPolynomial<R> operator * (const R c, const LaurentPolynomial<R>& p)
  {
    return (LaurentPolynomial<R>(p) *= c);
  }

//   template <class R> LaurentPolynomial<R> operator * (const LaurentPolynomial<R>& p, const LaurentPolynomial<R>& q)
//   {
//     LaurentPolynomial<R> r(p);
//     r *= q;

//     return r;
//   }

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

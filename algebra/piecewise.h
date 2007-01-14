// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_PIECEWISE_H
#define _MATHTL_PIECEWISE_H

#include <iostream>
#include <utils/function.h>
#include <algebra/polynomial.h>

namespace MathTL
{
  /*!
    A template class for univariate compactly supported functions f which
    admit a piecewise polynomial representation.

    The nodes of a piecewise function are assumed to take dyadic values
    2^{-j}k, with
    - supp(f) = [2^{-j}tmin, 2^{-j}tmax]
    - f|[2^{-j}k,2^{-j}(k+1)] is a polynomial for each integer k

    Note:
    - the local polynomial degree may vary
    - f may not be as regular as a spline

    Piecewise polynomial functions may be
    - constructed as B-splines
    - evaluated
    - added/subtracted
    - multiplied with constants
    - multiplied with polynomials
    - multiplied with another piecewise function
    - (symbolically) differentiated
    - (symbolically) integrated, quadrature rules optional

    We derive Piecewise<C> from the class Function<1>, as it is indeed one.
  */
  template <class C>
  class Piecewise
    : public Function<1, C>
  {
    public:
    typedef typename std::map<int, Polynomial<C> > PiecesType;

    /* data fields *********************************************************/
    protected:
    PiecesType expansion;  // the local polynomials
    int granularity;       // granularity j

    /* constructors and destructors ****************************************/
    public:
    /*!
      default constructor: yields the zero spline, granularity j=0
    */
    Piecewise();

    /*!
      copy constructor
    */
    Piecewise(const Piecewise<C>& p);

    /*!
      constructor with predefined granularity j
    */
    Piecewise(const int j);

    /*!
      destructor
    */
    ~Piecewise();

    /* member functions ****************************************************/
    /*!
      read access, local expansion
      get the local polynomial representation at interval k
    */
    Polynomial<C> get_local_expansion(const int k) const;

    /*!
      read access, all polynomials
      get all polynomial representations
    */
    const PiecesType* get_expansion() const;

    /*!
      write access, local expansion
      set the local polynomial representation at interval k
    */
    void set_local_expansion(const int k, const Polynomial<C> &p);

    /*!
      get the granularity j
    */
    inline int get_granularity() const { return granularity; }

    /*!
      clip this spline to [2^{-j}k1,2^{-j}k2]
    */
    void clip_me(const int k1, const int k2);

    /*!
      clip to [2^{-j}k1,2^{-j}k2]
    */
    Piecewise<C> clip(const int k1, const int k2) const;

    /*!
      increase granularity, jnew >= current granularity
    */
    void split_me(const int jnew);

    /*!
      increase granularity, jnew >= current granularity
    */  
    Piecewise split(const int jnew) const;

    /*!
      dilate by 2^{-j}
    */
    void dilate_me(const int j);

    /*!
      dilate by 2^{-j}
      (makes a copy)
    */
    Piecewise<C> dilate(const int j) const;

    /*!
      shift by 2^{-j}k
    */
    void shift_me(const int k);

    /*!
      shift by 2^{-j}k
      (makes a copy)
    */
    Piecewise<C> shift(const int k) const;

    /*!
      symbolic differentiation
    */
    Piecewise<C> differentiate() const;

    /*!
      integration over entire support
    */
    double integrate(const bool quadrature = false) const;

    /*!
      integration over [2^{-j}k1,2^{-j}k2]
    */
    double integrate(const int k1, const int k2,
                     const bool quadrature = false) const;

    /*!
      point evaluation
    */
    C value (const C x) const;

    /*!
      point evaluation
      (calls the above value(const C))
      needed as inheritance from Function<C>
    */
    inline C value (const Point<1>& x, const unsigned int component = 0) const;

    /*!
      point evaluation
      (calls the above value(const C))
      needed as inheritance from Function<C>
    */
    void vector_value(const Point<1> &p,
                      Vector<C>& values) const;

    /*!
      point evaluation of first derivative
    */
    C derivative (const C x) const;

    /*!
      in-place multiplication with a constant
    */
    Piecewise<C>& scale (const C c);

    /*!
      add a polynomial to this piecewise
    */
    Piecewise<C>& add (const Polynomial<C> &p);

    /*!
      add an other piecewise to this piecewise
    */
    Piecewise<C>& add (const Piecewise<C> &p);

    /*!
      inner product with another piecewise
      (don't use this extensively, since one copy has to be made!)
    */
    C inner_product(const Piecewise<C> &p) const;

    /*!
      inner product with another polynomial
      (don't use this extensively, since one copy has to be made!)
    */
    C inner_product(const Polynomial<C> &p) const;


    /* operators ***********************************************************/
    /*!
      assignment of another piecewise
    */
    Piecewise<C>& operator = (const Piecewise<C>& p);

    /*!
      point evaluation
    */
    C operator () (const C x) const;

    /*!
      in-place addition of a polynomial p
    */
    Piecewise<C>& operator += (const Polynomial<C> &p);

    /*!
      in-place addition of a piecewise p
    */
    Piecewise<C>& operator += (const Piecewise<C> &p);

    /*!
      subtraction of a polynomial
    */
    Piecewise<C>& operator -= (const Polynomial<C> &p);

    /*!
      subtraction of another piecewise
    */
    Piecewise<C>& operator -= (const Piecewise<C> &p);

    /*!
      scaling with a constant
    */
    Piecewise<C>& operator *= (const C c);

    /*!
      pointwise multiplication with a polynomial
    */
    Piecewise<C>& operator *= (const Polynomial<C> &p);

    /*!
      pointwise multiplication with another piecewise
    */
    Piecewise<C>& operator *= (const Piecewise<C> &p);

    /*!
      Error output
    */
    void MatError (char* str) const;
  };


  /* non class-member operators ********************************************/

  /*!
    pointwise sum of two polynomials
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  + (const Polynomial<C> &p, const Piecewise<C> &q);

  /*!
    pointwise sum of two polynomials
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  + (const Piecewise<C> &p, const Polynomial<C> &q);

  /*!
    pointwise sum of two polynomials
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  + (const Piecewise<C> &p, const Piecewise<C> &q);

  /*!
    pointwise difference of two piecewise
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  - (const Polynomial<C> &p, const Piecewise<C> &q);

  /*!
    pointwise difference of a piecewise from a polynomial
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  - (const Piecewise<C> &p, const Polynomial<C> &q);

  /*!
    pointwise difference of a polynomial from a piecewise
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  - (const Piecewise<C> &p, const Piecewise<C> &q);

  /*!
    multiplication with a real number
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  * (const C c, const Piecewise<C> &q);

  /*!
    pointwise multiplication with a polynomial
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  * (const Polynomial<C> &p, const Piecewise<C> &q);

  /*!
    pointwise multiplication with another piecewise
    (makes a copy)
  */
  template <class C>
  Piecewise<C> operator  * (const Piecewise<C> &p, const Piecewise<C> &q);

  /*!
    output into a stream
  */
  template <class C>
  std::ostream& operator << (std::ostream &s, const Piecewise<C>& p);
}

// implementation of inline functions
#include <algebra/piecewise.cpp>

#endif // _MATHTL_PIECEWISE_H

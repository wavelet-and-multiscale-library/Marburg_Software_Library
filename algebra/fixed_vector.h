// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_FIXED_VECTOR_H
#define _MATHTL_FIXED_VECTOR_H

#include <iostream>

// external functionality, for convenience:
#include <algebra/vector_norms.h>
#include <algebra/vector_arithmetics.h>

namespace MathTL
{
  /*!
    This class models finite, densely populated vectors
      x = (x_0, ... ,x_{n-1})
    of a a-priori known size with entries from an arbitrary (scalar)
    class C, designed for numerical computations.
    The signature parallels that of Vector.

    A FixedVector<C, SIZE> has essentially the same functionality as a
    FixedArray1D<C, SIZE>. However, we deliberately don't use inheritance
    here in order to really provide a raw class with maximum performance.
  */
  template <class C, unsigned int SIZE>
  class FixedVector
  {
  public:
    /*!
      value type (cf. STL containers)
    */
    typedef C value_type;
    
    /*!
      pointer type (cf. STL containers)
     */
    typedef value_type* pointer;

    /*!
      const pointer type (cf. STL containers)
    */
    typedef const value_type* const_pointer;

    /*!
      iterator type (cf. STL containers)
    */
    typedef value_type* iterator;

    /*!
      const iterator type (cf. STL containers)
    */
    typedef const value_type* const_iterator;

    /*!
      reference type (cf. STL containers)
    */
    typedef value_type& reference;
    
    /*!
      const reference type (cf. STL containers)
    */
    typedef const value_type& const_reference;

    /*!
      type of indexes and size type (cf. STL containers)
     */
    typedef unsigned int size_type;

    /*!
      construct a vector with all entries of the given number
      also serves as default constructor: yields a zero vector
     */
    FixedVector(const C c = 0);

    /*!
      copy constructor
    */
    FixedVector(const FixedVector<C, SIZE>& v);

    /*!
      release allocated memory
    */
    virtual ~FixedVector();

    /*!
      size/dimension of the vector
    */
    inline size_type size() const {return SIZE;}

    /*!
      assignment from another fixed vector
    */
    FixedVector<C, SIZE>& operator = (const FixedVector<C, SIZE>& v);

    /*!
      assignment of a constant value to each component
     */
    FixedVector<C, SIZE>& operator = (const C c);

    /*!
      C style read-only access to the i-th vector component
    */
    const C operator [] (const size_type i) const;

    /*!
      Matlab style read-only access to the i-th vector component
    */
    const C operator () (const size_type i) const;

    /*!
      C style read-write access to the i-th vector component
    */
    C& operator [] (const size_type i);

    /*!
      Matlab style read-write access to the i-th vector component
    */
    C& operator () (const size_type i);

    /*!
      read-only iterator access to first element (cf. STL containers)
    */
    const_iterator begin() const;

    /*!
      read-write iterator access to first element (cf. STL containers)
    */
    iterator begin();

    /*!
      read-only iterator access to the element behind the last one
      (cf. STL containers)
    */
    const_iterator end() const;

    /*!
      read-only iterator access to the element behind the last one
      (cf. STL containers)
    */
    iterator end();

    /*!
      equality test with another vector
     */
    template <class C2>
    bool operator == (const FixedVector<C2, SIZE>& v) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const FixedVector<C2, SIZE>& v) const;

    /*!
      lexicographical order
    */
    template <class C2>
    bool operator < (const FixedVector<C2, SIZE>& v) const;

    /*!
      in place summation *this += v
    */
    template <class C2>
    void add(const FixedVector<C2, SIZE>& v);
    
    /*!
      in place summation *this += s*v
    */
    template <class C2>
    void add(const C2 s, const FixedVector<C2, SIZE>& v);
    
    /*!
      in place summation *this = s*(*this) + v
      (AXPY level 1 BLAS routine)
    */
    template <class C2>
    void sadd(const C s, const FixedVector<C2, SIZE>& v);

    /*!
      in place scaling *this *= s
    */
    void scale(const C s);

    /*!
      in place summation
    */
    template <class C2>
    FixedVector<C, SIZE>& operator += (const FixedVector<C2, SIZE>& v);

    /*!
      in place subtraction *this -= v
    */
    template <class C2>
    void subtract(const FixedVector<C2, SIZE>& v);

    /*!
      in place subtraction
    */
    template <class C2>
    FixedVector<C, SIZE>& operator -= (const FixedVector<C2, SIZE>& v);

    /*!
      inner product
    */
    template <class C2>
    const C inner_product (const FixedVector<C2, SIZE>& v) const;

    /*!
      in place multiplication with a scalar
    */
    FixedVector<C, SIZE>& operator *= (const C c);
    
    /*!
      in place division by a (nontrivial) scalar
    */
    FixedVector<C, SIZE>& operator /= (const C c);
    
  protected:
    /*!
      internal storage is just a pointer to a C array
    */
    C* values_;
  };

  /*!
    sum of two vectors
    (you should avoid using this operator, since it requires one vector
    to be copied. Use += or add() instead!)
   */
  template <class C, class C2, unsigned int SIZE>
  FixedVector<C, SIZE> operator + (const FixedVector<C, SIZE>& v1, const FixedVector<C2, SIZE>& v2);

  /*!
    difference of two vectors
    (you should avoid using this operator, since it requires one vector
    to be copied. Use -= or sadd() instead!)
   */
  template <class C, class C2, unsigned int SIZE>
  FixedVector<C, SIZE> operator - (const FixedVector<C, SIZE>& v1, const FixedVector<C2, SIZE>& v2);

  /*!
    stream output for vectors
    (we overload this for code clarity)
  */
  template <class C, unsigned int SIZE>
  std::ostream& operator << (std::ostream& os, const FixedVector<C, SIZE>& v);
}

// include implementation of inline functions
#include <algebra/fixed_vector.cpp>

#endif // _MATHTL_FIXED_VECTOR_H

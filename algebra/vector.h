// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_VECTOR_H
#define _MATHTL_VECTOR_H

#include <iostream>

// external functionality, for convenience:
#include <utils/vector_norms.h>
#include <utils/vector_arithmetics.h>

namespace MathTL
{
  /*!
    This class models finite, densely populated vectors
      x = (x_0, ... ,x_{n-1})
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
    The signature parallels that of std::vector.

    A Vector<C> has essentially the same functionality as an Array1D<C>.
    However, we deliberately don't use inheritance here in order to
    really provide a raw class with maximum performance.
  */
  template <class C>
  class Vector
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
    typedef size_t size_type;

    /*!
      default constructor: yields a 0-dimensional vector
     */
    Vector();

    /*!
      Construct a vector of size s.
      It is possible to choose whether the vector will be initialized
      by zeros or not. For the initialization we use a constructor
      of the form C(0), which is available for the builtin types.
      The default behaviour is to indeed initialize the vector.
    */
    explicit Vector(const size_type s, const bool initialize = true);

    /*!
      copy constructor
    */
    Vector(const Vector<C>& v);

    /*!
      release allocated memory
    */
    ~Vector();

    /*!
      size/dimension of the vector
    */
    const size_type size() const;

    /*!
      resize vector, initialize with C(0) if desired
    */
    void resize(const size_type s, const bool initialize = true);

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
      assignment of a constant value to each component
     */
    Vector<C>& operator = (const C c);

    /*!
      assignment from another vector
    */
    Vector<C>& operator = (const Vector<C>& v);

    /*!
      swap components of two vectors
    */
    void swap (Vector<C>& v);

    /*!
      equality test with another vector
     */
    template <class C2>
    bool operator == (const Vector<C2>& v) const;

    /*!
      non-equality test
    */
    template <class C2>
    bool operator != (const Vector<C2>& v) const;

    /*!
      lexicographical order
    */
    template <class C2>
    bool operator < (const Vector<C2>& v) const;

    /*!
      in place summation *this += v
    */
    template <class C2>
    void add(const Vector<C2>& v);
    
    /*!
      in place summation *this = s*(*this) + v
      (AXPY level 1 BLAS routine)
    */
    template <class C2>
    void sadd(const C s, const Vector<C2>& v);

    /*!
      in place scaling *this *= s
    */
    void scale(const C s);

    /*!
      in place summation
    */
    template <class C2>
    Vector<C>& operator += (const Vector<C2>& v);

    /*!
      in place subtraction
    */
    template <class C2>
    Vector<C>& operator -= (const Vector<C2>& v);

    /*!
      in place multiplication with a scalar
    */
    Vector<C>& operator *= (const C c);
    
    /*!
      in place division by a (nontrivial) scalar
    */
    Vector<C>& operator /= (const C c);

    /*!
      inner product
    */
    template <class C2>
    const C operator * (const Vector<C2>& v) const;

  protected:
    /*!
      internal storage is just a pointer to a C array
    */
    C* values_;

    /*!
      size/dimension of the vector
    */
    size_type size_;
  };

  /*!
    swap the values of two vectors
  */
  template <class C>
  void swap(Vector<C>& v1, Vector<C>& v2);
  
  /*!
    stream output for vectors
    (we overload this for code clarity)
  */
  template <class C>
  std::ostream& operator << (std::ostream& os, const Vector<C>& v);
}

// include implementation of inline functions
#include <algebra/vector.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_ARRAY1D_H
#define _MATHTL_ARRAY1D_H

#include <iostream>

namespace MathTL
{
  /*!
    This class models one-dimensional arrays of objects from an arbitrary
    class C and is merely a wrapper class around a conventional C-style array.
    Array1D<C> has the same interface as std::vector<C>, but it is much
    faster in practice due to the internal storage format.
    Even compared with a std::valarray<C>, one gains a bit of performance.
  */
  template <class C>
    class Array1D
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
      type of indexes and size of the array
     */
    typedef size_t size_type;

    /*!
      default constructor, yields an empty array
     */
    Array1D();

    /*!
      copy constructor
    */
    Array1D(const Array1D<C>& a);
    
    /*!
      Construct an array of positive size,
      please note that the builtin types (int, double, ...) may _not_
      automatically be set to _zero_, this depends on the compiler!
    */
    explicit Array1D(const size_type s);

    /*!
      release allocated memory
    */
    ~Array1D();

    /*!
      size of the array
    */
    const size_type size() const;

    /*!
      resize the array

      Warning: the initialization of the entries depends on the compiler (and even on the compiler options).
      It may very well be something like NAN!
    */
    void resize(const size_type s);

    /*!
      assignment
    */
    Array1D<C>& operator = (const Array1D<C>& a);

    /*!
      read-only access to the i-th array member
    */
    const C& operator [] (const size_type i) const;

    /*!
      read-write access to the i-th array member
    */
    C& operator [] (const size_type i);

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
      swap components of two arrays
    */
    void swap (Array1D<C>& a);

    /*!
      swap two entries of an array
    */
    void swap (const size_type i, const size_type j);

  protected:
    /*!
      internal storage is just a pointer to a C array
    */
    C* data_;

    /*!
      size of the array
    */
    size_type size_;
  };

  /*!
    stream output for arrays
   */
  template <class C>
  std::ostream& operator << (std::ostream& os, const Array1D<C>& A);
}

// include implementation of inline functions
#include "utils/array1d.cpp"

#endif

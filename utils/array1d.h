// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
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
      type of the objects stored in this array
     */
    typedef C value_type;

    /*!
      type of indexes and size of the array
     */
    typedef size_t size_type;

    /*!
      default constructor, yields an empty array
     */
    Array1D();

    /*!
      construct an array of positive size and initialize its elements
      (Please note that the builtin types (int, double, ...) are not
      automatically set to zero, this depends on the compiler!)
    */
    Array1D(const size_type s);

    /*!
      release allocated memory
    */
    ~Array1D();

    /*!
      size of the array
    */
    const size_type size() const;

    /*!
      read-only access to the i-th array member
    */
    const C& operator [] (const size_type i) const;

    /*!
      read-write access to the i-th array member
    */
    C& operator [] (const size_type i);

  private:
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

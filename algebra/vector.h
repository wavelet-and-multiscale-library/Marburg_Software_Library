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
#include "utils/array1d.h"

namespace MathTL
{
  /*!
    This class models finite, densely populated vectors
      x = (x_0, ... ,x_{n-1})
    with entries from an arbitrary (scalar) class C,
    designed for numerical computations.
    The signature parallels that of std::vector.

    A Vector<C> is derived from Array1D<C>, so we also refer to
    the Array1D<C> documentation here.
  */
  template <class C>
  class Vector : public Array1D<C>
  {
  public:
    /*!
      value type (cf. STL containers)
     */
    typedef typename Array1D<C>::value_type value_type;

    /*!
      pointer type (cf. STL containers)
     */
    typedef typename Array1D<C>::pointer pointer;

    /*!
      const pointer type (cf. STL containers)
    */
    typedef typename Array1D<C>::const_pointer const_pointer;

    /*!
      iterator type (cf. STL containers)
    */
    typedef typename Array1D<C>::iterator iterator;

    /*!
      const iterator type (cf. STL containers)
    */
    typedef typename Array1D<C>::const_iterator const_iterator;

    /*!
      reference type (cf. STL containers)
    */
    typedef typename Array1D<C>::reference reference;
    
    /*!
      const reference type (cf. STL containers)
    */
    typedef typename Array1D<C>::const_reference const_reference;

    /*!
      type of indexes and size type (cf. STL containers)
    */
    typedef typename Array1D<C>::size_type size_type;

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
      assignment of a constant value to each component
     */
    Vector<C>& operator = (const C& c);
  };

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

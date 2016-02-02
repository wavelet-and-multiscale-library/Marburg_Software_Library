// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_POINT_H
#define _MATHTL_POINT_H

#include <iostream>
#include "algebra/tensor.h"

namespace MathTL
{
  /*!
    a class for points in the d-dimensional Euclidean space
    where d is a priori known
  */
  template <unsigned int DIM, class VALUE = double>
  class Point : public Tensor<1, DIM, VALUE>
  {
  public:
    /*!
      value type of the vector (cf. STL containers)
     */
    typedef VALUE value_type;

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
      size type (cf. STL containers)
    */
    typedef size_t size_type;

    /*!
      default constructor: yields the origin
      (we always initialize the point with zero)
    */
    Point();
    
    /*!
      copy constructor from a tensor of rank 1
    */
    Point(const Tensor<1, DIM, VALUE>&);

    /*!
      constructor from a single real number, sets all coordinates to this number
    */
    explicit Point(const VALUE x);

    /*!
      constructor from two real numbers, this is only allowed for DIM==2
    */
    Point(const VALUE x, const VALUE y);

    /*!
      constructor from three real numbers, this is only allowed for DIM==3
    */
    Point(const VALUE x, const VALUE y, const VALUE z);

    /*!
      set all coordinates to a real number
      (also for implicit conversion)
    */
    Point<DIM, VALUE>& operator = (const VALUE x);

//     /*!
//       assignment operator (for safety, Tensor class already has one)
//     */
//     Point<DIM, VALUE>& operator = (const Point<DIM, VALUE>& p);

    /*!
      size/dimension of the point (cf. std::vector signature)
     */
    const size_type size() const;

    /*!
      Matlab-style read-only access operator (operator [] is inherited
      from the Tensor class)
    */
    const VALUE operator () (const size_type i) const;

    /*!
      Matlab-style read-write access operator (operator [] is inherited
      from the Tensor class)
    */
    VALUE& operator () (const size_type i);
  };

  /*!
    stream output for points (overloaded from the tensor one)
   */
  template <unsigned int DIM, class VALUE>
  std::ostream& operator << (std::ostream& os, const Point<DIM, VALUE>& p);
}

#include "geometry/point.cpp"

#endif

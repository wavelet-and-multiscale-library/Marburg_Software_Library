// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_POINT_H
#define _MATHTL_POINT_H

#include "algebra/tensor.h"

namespace MathTL
{
  // a class for points in the d-dimensional Euclidean space
  template <int DIM>
    class Point : public Tensor<1, DIM>
  {
  public:
    /*!
      value type of the vector (cf. STL containers)
     */
    typedef double value_type;

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
    Point(const Tensor<1, DIM>&);

    /*!
      constructor from a real number, this is only allowed for DIM==1
    */
    explicit Point(const double x);

    /*!
      constructor from two real numbers, this is only allowed for DIM==2
    */
    Point(const double x, const double y);

    /*!
      constructor from three real numbers, this is only allowed for DIM==3
    */
    Point(const double x, const double y, const double z);

    /*!
      Matlab-style read-only access operator (operator [] is inherited
      from the Tensor class)
    */
    const double operator () (const size_type i) const;

    /*!
      Matlab-style read-write access operator (operator [] is inherited
      from the Tensor class)
    */
    double& operator () (const size_type i);
  };
}

#include "geometry/point.cpp"

#endif

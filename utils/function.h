// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_FUNCTION_H
#define _MATHTL_FUNCTION_H

#include <algebra/vector.h>
#include <utils/function_time.h>
#include <geometry/point.h>

namespace MathTL
{
  template <unsigned int DIM, class VALUE> class Point;
  template <unsigned int RANK, unsigned int DIM, class VALUE> class Tensor;
  template <unsigned int DIM, class VALUE> class Tensor<1, DIM, VALUE>;

  /*!
    Base class for scalar or vector-valued functions on \mathbb R^d

    (this function model is essetially compatible to that of the deal.II
    library, version 5)
   */
  template <unsigned int DIM, class VALUE = double>
  class Function : public FunctionTime
  {
  public:
    /*!
      make template parameter value accessible
     */
    static const unsigned int dimension = DIM;

    /*!
      size type
    */
    typedef size_t size_type;

    /*!
      dimension of value
    */
    const unsigned int n_components;

    /*!
      default constructor:
      with the default values, this yields a scalar function
      with initial time value 0.
     */
    Function(const unsigned int n_components = 1,
	     const double initial_time = 0.0);

    /*!
      purely virtual destructor
    */
    virtual ~Function() = 0;

    /*!
      evaluate (a component of) the function
    */
    virtual VALUE value (const Point<DIM>& p,
			 const unsigned int component = 0) const = 0;

    /*!
      evaluate the function
    */
    virtual void vector_value(const Point<DIM> &p,
			      Vector<VALUE>& values) const = 0;

    /*!
      (estimate of) memory consumption in bytes
    */
    const size_type memory_consumption() const;
  };


  //
  //
  // some example functions

  /*!
    scalar or vector valued zero function
   */
  template <unsigned int DIM, class VALUE = double>
  class ZeroFunction
    : public Function<DIM, VALUE>
  {
  public:
    /*!
      default constructor, one component (scalar-valued) by default
    */
    ZeroFunction(const unsigned int n_components = 1);

    /*!
      virtual (!) destructor
    */
    virtual ~ZeroFunction();

    /*!
      evaluate one component of the zero function
     */
    VALUE value(const Point<DIM>& p,
		const unsigned int component = 0) const;

    /*!
      evaluate the zero function
    */
    void vector_value(const Point<DIM> &p,
		      Vector<VALUE>& values) const;
  };
  
}

// implementations of inline functions
#include "utils/function.cpp"

#endif

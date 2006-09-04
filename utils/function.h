// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
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

    (this function model is essentially compatible to that of the deal.II
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
    virtual VALUE value (const Point<DIM,VALUE>& p,
			 const unsigned int component = 0) const = 0;

    /*!
      evaluate the function
      (values should be of appropriate size)
    */
    virtual void vector_value(const Point<DIM,VALUE> &p,
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
    ZeroFunction(const unsigned int n_components = 1);
    virtual ~ZeroFunction();
    VALUE value(const Point<DIM,VALUE>& p,
		const unsigned int component = 0) const;
    void vector_value(const Point<DIM,VALUE> &p,
		      Vector<VALUE>& values) const;
  };

  /*!
    constant function
  */
  template <unsigned int DIM, class VALUE = double>
  class ConstantFunction
    : public Function<DIM, VALUE>
  {
  public:
    ConstantFunction(const Vector<VALUE>& value);
    virtual ~ConstantFunction();
    VALUE value(const Point<DIM,VALUE>& p,
		const unsigned int component = 0) const;
    void vector_value(const Point<DIM,VALUE> &p,
		      Vector<VALUE>& values) const;
  protected:
    //! the value
    Vector<VALUE> c;
  };

  /*!
    pointwise product of two (preferably scalar--valued) functions
  */
  template <unsigned int DIM, class VALUE = double>
  class ProductFunction
    : public Function<DIM, VALUE>
  {
  public:
    ProductFunction(const Function<DIM,VALUE>* f1,
		    const Function<DIM,VALUE>* f2);
    virtual ~ProductFunction();
    VALUE value(const Point<DIM,VALUE>& p,
		const unsigned int component = 0) const;
    void vector_value(const Point<DIM,VALUE> &p,
		      Vector<VALUE>& values) const;
    
  protected:
    const Function<DIM,VALUE>* f1_;
    const Function<DIM,VALUE>* f2_;
  };
}

// implementations of inline functions
#include "utils/function.cpp"

#endif

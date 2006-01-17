// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_TENSOR_BASE_H
#define _MATHTL_TENSOR_BASE_H

#include <iostream>
#include "algebra/tensor.h"
#include "utils/array1d.h"

namespace MathTL
{
  // general template declaration, see tensor.h
  template <unsigned int RANK, unsigned int DIM, class VALUE> class Tensor;
  template <unsigned int DIM, class VALUE> class Tensor<1, DIM, VALUE>;

  /*!
    Special version of tensors of rank 1 on the Euclidean space
    \mathbb R^d, i.e., vectors of dimension DIM.
  */
  template <unsigned int DIM, class VALUE>
  class Tensor<1, DIM, VALUE>
  {
  public:
    /*!
      type of internal storage object: a real number
     */
    typedef VALUE value_type;

    /*!
      size type
    */
    typedef size_t size_type;

    /*!
      dimension of the tensor, available at compile time
     */
    static const unsigned int dimension = DIM;

    /*!
      rank of the tensor, available at compile time
    */
    static const unsigned int rank = 1;

    /*!
      default constructor, yields a zero tensor;
      it is possible to choose whether the tensor elements are
      initialized by zero (default behaviour) or not
    */
    explicit Tensor(const bool initialize = true);

    /*!
      reset Tensor to a zero one,
      this is equivalent to the state given by the default constructor
    */
    void clear();

    /*!
      size of the tensor as a vector (cf. std::vector signature)
    */
    const size_type size() const;
    
    /*!
      read-only access to the i-th coordinate
    */
    const VALUE operator [] (const size_type i) const;

    /*!
      read-write access to the i-th coordinate
    */
    VALUE& operator [] (const size_type i);

    /*!
      equality test of two tensors
    */
    bool operator == (const Tensor<1, DIM, VALUE>&) const;

    /*!
      non-equality test of two tensors
    */
    bool operator != (const Tensor<1, DIM, VALUE>&) const;

    /*!
      in-place summation of two tensors
    */
    Tensor<1, DIM, VALUE>& operator += (const Tensor<1, DIM, VALUE>&);

    /*!
      in place summation *this += s*t
    */
    void add(const VALUE s, const Tensor<1, DIM, VALUE>& T);

    /*!
      in-place subtraction of two tensors
    */
    Tensor<1, DIM, VALUE>& operator -= (const Tensor<1, DIM, VALUE>&);

    /*!
      in-place multiplication of a tensor with a scalar from the left
    */
    Tensor<1, DIM, VALUE>& operator *= (const VALUE s);

    /*!
      in-place division of a tensor by a scalar
    */
    Tensor<1, DIM, VALUE>& operator /= (const VALUE s);

    /*!
      sum of two tensors (makes a copy)
    */
    Tensor<1, DIM, VALUE> operator + (const Tensor<1, DIM, VALUE>&) const;

    /*!
      difference of two tensors (makes a copy)
    */
    Tensor<1, DIM, VALUE> operator - (const Tensor<1, DIM, VALUE>&) const;

    /*!
      tensor with negated entries
    */
    Tensor<1, DIM, VALUE> operator - () const;
    
    /*!
      inner product of two rank 1 tensors (scalar product)
    */
    VALUE operator * (const Tensor<1, DIM, VALUE>&) const;

    /*!
      estimate memory consumption in bytes
    */
    static const size_type memory_consumption();
    
  protected:
    /*!
      internal storage for the vector entries
    */
    Array1D<VALUE> values;
  };

  /*!
    stream output for rank 1 tensors
  */
  template <unsigned int DIM, class VALUE>
    std::ostream& operator << (std::ostream& os,
			       const Tensor<1, DIM, VALUE>& T);

  /*!
    stream output for rank and dimension 1 tensors
  */
  template <class VALUE>
  std::ostream& operator << (std::ostream& os,
			     const Tensor<1, 1, VALUE>& T);
}

#include <algebra/tensor_base.cpp>

#endif

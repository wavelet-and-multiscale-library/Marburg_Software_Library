// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_TENSOR_BASE_H
#define _MATHTL_TENSOR_BASE_H

#include <iostream>
#include "algebra/tensor.h"
#include "utils/array1d.h"

namespace MathTL
{
  // general template declaration, see tensor.h
  template <unsigned int RANK, unsigned int DIM> class Tensor;
  template <unsigned int DIM> class Tensor<1, DIM>;

  /*!
    Special version of tensors of rank 1 on the Euclidean space
    \mathbb R^d, i.e., vectors of dimension DIM.
  */
  template <unsigned int DIM>
    class Tensor<1, DIM>
  {
  public:
    /*!
      type of internal storage object: a real number
     */
    typedef double value_type;

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
      read-only access to the i-th coordinate
    */
    const double operator [] (const size_type i) const;

    /*!
      read-write access to the i-th coordinate
    */
    double& operator [] (const size_type i);

    /*!
      equality test of two tensors
    */
    bool operator == (const Tensor<1, DIM>&) const;

    /*!
      non-equality test of two tensors
    */
    bool operator != (const Tensor<1, DIM>&) const;

    /*!
      in-place addition of two tensors
    */
    Tensor<1, DIM>& operator += (const Tensor<1, DIM>&);

    /*!
      in-place subtraction of two tensors
    */
    Tensor<1, DIM>& operator -= (const Tensor<1, DIM>&);

    /*!
      in-place multiplication of a tensor with a scalar from the left
    */
    Tensor<1, DIM>& operator *= (const double s);

    /*!
      in-place division of a tensor by a scalar
    */
    Tensor<1, DIM>& operator /= (const double s);

    /*!
      sum of two tensors (makes a copy)
    */
    Tensor<1, DIM> operator + (const Tensor<1, DIM>&) const;

    /*!
      difference of two tensors (makes a copy)
    */
    Tensor<1, DIM> operator - (const Tensor<1, DIM>&) const;

    /*!
      tensor with negated entries
    */
    Tensor<1, DIM> operator - () const;
    
    /*!
      inner product of two rank 1 tensors (scalar product)
    */
    double operator * (const Tensor<1, DIM>&) const;

    /*!
      estimate memory consumption in bytes
    */
    static const size_type memory_consumption();
    
  protected:
    /*!
      storage for the vector entries
      (TODO: use my own array class)
    */
    Array1D<double> values;
  };

  /*!
    stream output for rank 1 tensors
  */
  template <unsigned int DIM>
    std::ostream& operator << (std::ostream& os, const Tensor<1, DIM>& T);

  /*!
    stream output for rank and dimension 1 tensors
  */
  std::ostream& operator << (std::ostream& os, const Tensor<1, 1>& T);

}

#include <algebra/tensor_base.cpp>

#endif

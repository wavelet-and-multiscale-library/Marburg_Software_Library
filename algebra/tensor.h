// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_TENSOR_H
#define _MATHTL_TENSOR_H

#include <iostream>
#include "algebra/tensor_base.h"

namespace MathTL
{
  template <unsigned int RANK, unsigned int DIM, class VALUE> class Tensor;
  template <unsigned int DIM, class VALUE> class Tensor<1, DIM, VALUE>;

  /*!
    This class models tensors of an arbitrary rank r on the Euclidean space
    \mathbb R^d. Most functionality can be handed over to tensors of lower-dimensional
    rank. Tensors of rank 1 are modeled via template specialization
    in tensor_base.h .
  */
  template <unsigned int RANK, unsigned int DIM, class VALUE = double>
  class Tensor
  {
  public:
    /*!
      type of internal storage object, due to the recursive structure it
      is just a tensor with lower rank
     */
    typedef Tensor<RANK-1, DIM, VALUE> value_type;

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
    static const unsigned int rank = RANK;

    /*!
      default constructor, yields a zero tensor
    */
    Tensor();

    /*!
      reset Tensor to a zero one,
      this is equivalent to the state given by the default constructor
    */
    void clear();

    /*!
      read-only access to the subtensors
    */
    const Tensor<RANK-1, DIM, VALUE>& operator [] (const size_type i) const;

    /*!
      read-write access to the subtensors
    */
    Tensor<RANK-1, DIM, VALUE>& operator [] (const size_type i);

    /*!
      assignment of a tensor
    */
    Tensor<RANK, DIM, VALUE>& operator = (const Tensor<RANK, DIM, VALUE>&);

    /*!
      equality test of two tensors
    */
    bool operator == (const Tensor<RANK, DIM, VALUE>&) const;

    /*!
      non-equality test of two tensors
    */
    bool operator != (const Tensor<RANK, DIM, VALUE>&) const;

    /*!
      in-place sum of two tensors
    */
    Tensor<RANK, DIM, VALUE>& operator += (const Tensor<RANK, DIM, VALUE>&);

    /*!
      in place summation *this += s*t
    */
    void add(const VALUE s, const Tensor<RANK, DIM, VALUE>& T);
    
    /*!
      in-place difference of two tensors
    */
    Tensor<RANK, DIM, VALUE>& operator -= (const Tensor<RANK, DIM, VALUE>&);

    /*!
      in-place multiplication of a tensor with a scalar from the left
    */
    Tensor<RANK, DIM, VALUE>& operator *= (const VALUE s);

    /*!
      in-place division of a tensor by a scalar
    */
    Tensor<RANK, DIM, VALUE>& operator /= (const VALUE s);

    /*!
      sum of two tensors (makes a copy)
    */
    Tensor<RANK, DIM, VALUE> operator + (const Tensor<RANK, DIM, VALUE>&) const;

    /*!
      difference of two tensors (makes a copy)
    */
    Tensor<RANK, DIM, VALUE> operator - (const Tensor<RANK, DIM, VALUE>&) const;

    /*!
      tensor with negated entries
    */
    Tensor<RANK, DIM, VALUE> operator - () const;
    
    /*!
      estimate memory consumption in bytes
    */
    static const size_type memory_consumption();
    
  protected:
    // internally, we store a tensor of rank r just as a d-dimensional array
    // of tensors of lower rank
    Tensor<RANK-1, DIM, VALUE> subtensor[DIM];
  };

  //
  //
  // external functionality for general tensors

  /*!
    stream output for general tensors:
    due to the recursive structure, the elements of the tensor are printed
    with one blank in between, two blanks between rank 1 subtensors and so on.
  */
  template <unsigned int RANK, unsigned int DIM, class VALUE>
  std::ostream& operator << (std::ostream& os,
			     const Tensor<RANK, DIM, VALUE>& T);

  /*!
    stream output for dimension 1 tensors
  */
  template <unsigned int RANK, class VALUE>
    std::ostream& operator << (std::ostream& os,
			       const Tensor<RANK, 1, VALUE>& T);
}

// include implementation of inline functions
#include <algebra/tensor.cpp>

#endif

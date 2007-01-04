// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_FIXED_ARRAY1D_H
#define _MATHTL_FIXED_ARRAY1D_H

#include <iostream>

namespace MathTL
{
  /*!
    This class models one-dimensional arrays of objects from an arbitrary
    class C, where the array size is a priori known.
    The signature is identical to Array1D<C>.
  */
  template <class C, unsigned int SIZE>
  class FixedArray1D
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
    typedef unsigned int size_type;
    
    /*!
      default constructor, yields an empty array
    */
    FixedArray1D();
    
    /*!
      copy constructor
    */
    FixedArray1D(const FixedArray1D<C, SIZE>& a);
    
    /*!
      release allocated memory
    */
    ~FixedArray1D();
    
    /*!
      size of the array
    */
    const size_type size() const;
    
    /*!
      assignment
    */
    FixedArray1D<C, SIZE>& operator = (const FixedArray1D<C, SIZE>& a);
    
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
    
  protected:
    /*!
      internal storage is just a pointer to a C array
    */
    C* data_;
  };

  //! template specialization to C=double
  template <unsigned int SIZE>
  class FixedArray1D<double,SIZE>
  {
  public:
    //! value type (cf. STL containers)
    typedef double value_type;
    
    //! pointer type (cf. STL containers)
    typedef value_type* pointer;
    
    //! const pointer type (cf. STL containers)
    typedef const value_type* const_pointer;
    
    //! iterator type (cf. STL containers)
    typedef value_type* iterator;
    
    //! const iterator type (cf. STL containers)
    typedef const value_type* const_iterator;
    
    //! reference type (cf. STL containers)
    typedef value_type& reference;
    
    //! const reference type (cf. STL containers)
    typedef const value_type& const_reference;
    
    //! type of indexes and size of the array
    typedef unsigned int size_type;
    
    //! default constructor, yields an empty array
    FixedArray1D();
    
    //! copy constructor
    FixedArray1D(const FixedArray1D<double, SIZE>& a);
    
    //! release allocated memory
    ~FixedArray1D();
    
    //! size of the array
    const size_type size() const;
    
    //! assignment
    FixedArray1D<double, SIZE>& operator = (const FixedArray1D<double, SIZE>& a);
    
    //! read-only access to the i-th array member
    const double& operator [] (const size_type i) const;
    
    //! read-write access to the i-th array member
    double& operator [] (const size_type i);
    
    //! read-only iterator access to first element (cf. STL containers)
    const_iterator begin() const;
    
    //! read-write iterator access to first element (cf. STL containers)
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
    
  protected:
    /*!
      internal storage is just a pointer to a C array,
      the "+1" is a temporary hack to allow for 0-dim. arrays
    */
    double data_[SIZE+1];
  };

  
  /*!
    stream output for arrays
   */
  template <class C, unsigned int SIZE>
  std::ostream& operator << (std::ostream& os, const FixedArray1D<C, SIZE>& A);
}

// include implementation of inline functions
#include "utils/fixed_array1d.cpp"

#endif

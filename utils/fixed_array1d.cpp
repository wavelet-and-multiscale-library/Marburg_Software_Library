// implementation of some (inline) FixedArray1D<C, SIZE>:: methods

#include <cassert>
#include <algorithm>
#include "io/vector_io.h"

namespace MathTL
{
  template <class C, unsigned int SIZE>
  inline
  FixedArray1D<C, SIZE>::FixedArray1D()
  {
    data_ = new C[SIZE];
  }
  
  template <class C, unsigned int SIZE>
  inline
  FixedArray1D<C, SIZE>::FixedArray1D(const FixedArray1D<C, SIZE>& a)
  {
    if (SIZE > 0)
      {
	data_ = new C[SIZE];
	std::copy(a.begin(), a.end(), begin());
      }
  }
  
  template <class C, unsigned int SIZE>
  FixedArray1D<C, SIZE>&
  FixedArray1D<C, SIZE>::operator = (const FixedArray1D<C, SIZE>& a)
  {
    std::copy(a.begin(), a.end(), begin());
    
    return *this;
  }
  
  template <class C, unsigned int SIZE>
  inline
  FixedArray1D<C, SIZE>::~FixedArray1D()
  {
    delete [] data_;
  }
  
  template <class C, unsigned int SIZE>
  inline
  const typename FixedArray1D<C, SIZE>::size_type
  FixedArray1D<C, SIZE>::size() const
  {
    return SIZE;
  }

  template <class C, unsigned int SIZE>
  inline
  const C& FixedArray1D<C, SIZE>::operator [] (const size_type i) const
  {
    assert(i < SIZE);
    return data_[i];
  }

  template <class C, unsigned int SIZE>
  inline
  C& FixedArray1D<C, SIZE>::operator [] (const size_type i)
  {
    assert(i < SIZE);
    return data_[i];
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedArray1D<C, SIZE>::const_iterator
  FixedArray1D<C, SIZE>::begin() const
  {
    return &data_[0];
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedArray1D<C, SIZE>::iterator
  FixedArray1D<C, SIZE>::begin()
  {
    return &data_[0];
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedArray1D<C, SIZE>::const_iterator
  FixedArray1D<C, SIZE>::end() const
  {
    return &data_[SIZE];
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedArray1D<C, SIZE>::iterator
  FixedArray1D<C, SIZE>::end()
  {
    return &data_[SIZE];
  }

  template <class C, unsigned int SIZE>
  inline
  std::ostream& operator << (std::ostream& os, const FixedArray1D<C, SIZE>& A)
  {
    print_vector(A, os);
    return os;
  }
}

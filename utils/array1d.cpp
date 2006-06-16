// implementation of some (inline) Array1D<C>:: methods

#include <cassert>
#include <algorithm>
#include "io/vector_io.h"

namespace MathTL
{
  template <class C>
  inline
  Array1D<C>::Array1D()
    : data_(0), size_(0)
  {
  }

  template <class C>
  inline
  Array1D<C>::Array1D(const size_type s)
    : size_(s)
  {
    if (s == 0)
      data_ = 0;
    else
      {
	data_ = new C[s]; // calls C()
      }
  }

  template <class C>
  inline
  Array1D<C>::Array1D(const Array1D<C>& a)
    : size_(a.size())
  {
    if (size_ == 0)
      data_ = 0;
    else
      {
	data_ = new C[size_];
	std::copy(a.begin(), a.end(), begin());
      }
  }

  template <class C>
  Array1D<C>& Array1D<C>::operator = (const Array1D<C>& a)
  {
    resize(a.size());
    std::copy(a.begin(), a.end(), begin());

    return *this;
  }

  template <class C>
  inline
  Array1D<C>::~Array1D()
  {
    if (data_ != 0) delete [] data_;
    size_ = 0;
  }
  
  template <class C>
  inline
  const typename Array1D<C>::size_type
  Array1D<C>::size() const
  {
    return size_;
  }

  template <class C>
  void Array1D<C>::resize(const size_type s)
  {
    if (s == 0)
      {
	if (data_ != 0) {
	  delete [] data_;
	  data_ = 0;
	}
	size_ = 0;
      }
    else
      {
	if (size_ != s) {
	  if (data_ != 0) delete [] data_;
	  data_ = new C[s]; // calls C()
	  size_ = s;
	}
      }
  }

  template <class C>
  inline
  const C& Array1D<C>::operator [] (const size_type i) const
  {
    assert(i < size_);
    return data_[i];
  }

  template <class C>
  inline
  C& Array1D<C>::operator [] (const size_type i)
  {
    assert(i < size_);
    return data_[i];
  }

  template <class C>
  inline
  typename Array1D<C>::const_iterator
  Array1D<C>::begin() const
  {
    return &data_[0];
  }

  template <class C>
  inline
  typename Array1D<C>::iterator
  Array1D<C>::begin()
  {
    return &data_[0];
  }

  template <class C>
  inline
  typename Array1D<C>::const_iterator
  Array1D<C>::end() const
  {
    return &data_[size_];
  }

  template <class C>
  inline
  typename Array1D<C>::iterator
  Array1D<C>::end()
  {
    return &data_[size_];
  }

  template <class C>
  inline
  void Array1D<C>::swap(Array1D<C>& a)
  {
    std::swap(data_, a.data_);
  }

  template <class C>
  inline
  std::ostream& operator << (std::ostream& os, const Array1D<C>& A)
  {
    print_vector(A, os);
    return os;
  }
}

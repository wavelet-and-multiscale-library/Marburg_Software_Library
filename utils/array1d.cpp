// implementation of some (inline) Array1D<C>:: methods

#include <cassert>

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
      data_ = new C[s]; // calls C()
  }

  template <class C>
  inline
  Array1D<C>::~Array1D()
  {
    if (data_ != 0)
      delete [] data_;
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
}

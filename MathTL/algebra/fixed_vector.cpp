// implementation of some (inline) FixedVector<C, SIZE>:: methods

#include <cassert>
#include <algorithm>
#include "io/vector_io.h"

namespace MathTL
{

   template <class C, unsigned int SIZE>
   inline
   FixedVector<C, SIZE>::FixedVector()
   {
    values_ = new C[SIZE];
    (*this).operator = (C());
   }

   template <class C, unsigned int SIZE>
   inline
   FixedVector<C, SIZE>::FixedVector(const C c)
   {
    values_ = new C[SIZE];
    (*this).operator = (c);
   }
  
  template <class C, unsigned int SIZE>
  inline
  FixedVector<C, SIZE>::FixedVector(const FixedVector<C, SIZE>& a)
  {
    if (SIZE > 0)
      {
        values_ = new C[SIZE];
        std::copy(a.begin(), a.end(), begin());
      }
  }
  
  template <class C, unsigned int SIZE>
  inline
  FixedVector<C, SIZE>::~FixedVector()
  {
    delete [] values_;
  }
  
  template <class C, unsigned int SIZE>
  FixedVector<C, SIZE>&
  FixedVector<C, SIZE>::operator = (const FixedVector<C, SIZE>& a)
  {
    std::copy(a.begin(), a.end(), begin());
    
    return *this;
  }
  
  template <class C, unsigned int SIZE>
  inline
  FixedVector<C, SIZE>&
  FixedVector<C, SIZE>::operator = (const C c)
  {
    assert(size() > 0);

    std::fill(begin(), end(), c);

    return *this;
  }

  template <class C, unsigned int SIZE>
  inline
  const C FixedVector<C, SIZE>::operator [] (const size_type i) const
  {
    assert(i < SIZE);
    return values_[i];
  }

  template <class C, unsigned int SIZE>
  inline
  const C FixedVector<C, SIZE>::operator () (const size_type i) const
  {
    return operator[] (i); // handles the assertions
  }

  template <class C, unsigned int SIZE>
  inline
  C& FixedVector<C, SIZE>::operator [] (const size_type i)
  {
    assert(i < SIZE);
    return values_[i];
  }

  template <class C, unsigned int SIZE>
  inline
  C& FixedVector<C, SIZE>::operator () (const size_type i)
  {
    return operator[] (i); // handles the assertions
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedVector<C, SIZE>::const_iterator
  FixedVector<C, SIZE>::begin() const
  {
    return &values_[0];
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedVector<C, SIZE>::iterator
  FixedVector<C, SIZE>::begin()
  {
    return &values_[0];
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedVector<C, SIZE>::const_iterator
  FixedVector<C, SIZE>::end() const
  {
    return &values_[SIZE];
  }

  template <class C, unsigned int SIZE>
  inline
  typename FixedVector<C, SIZE>::iterator
  FixedVector<C, SIZE>::end()
  {
    return &values_[SIZE];
  }

  template <class C, unsigned int SIZE>
  template <class C2>
  bool FixedVector<C, SIZE>::operator == (const FixedVector<C2, SIZE>& v) const
  {
    if (SIZE != v.size()) return false;
    return std::equal(begin(), end(), v.begin());
  }

  template <class C, unsigned int SIZE>
  template <class C2>
  inline
  bool FixedVector<C, SIZE>::operator != (const FixedVector<C2, SIZE>& v) const
  {
    return !((*this) == v);
  }

  template <class C, unsigned int SIZE>
  template <class C2>
  inline
  bool FixedVector<C, SIZE>::operator < (const FixedVector<C2, SIZE>& v) const
  {
    assert(SIZE == v.size());

    return std::lexicographical_compare(begin(), end(), v.begin(), v.end());
  }

  template <class C, unsigned int SIZE>
  template <class C2>
  void FixedVector<C, SIZE>::add(const FixedVector<C2, SIZE>& v)
  {
    assert(SIZE > 0);
    assert(SIZE == v.size());

    iterator it(begin()), itend(end());
    typename FixedVector<C2, SIZE>::const_iterator itv(v.begin());
    while (it != itend)
      *it++ += *itv++;
  }
   
  template <class C, unsigned int SIZE>
  template <class C2>
  void FixedVector<C, SIZE>::add(const C2 s, const FixedVector<C2, SIZE>& v)
  {
    assert(SIZE > 0);
    assert(SIZE == v.size());

    iterator it(begin()), itend(end());
    typename FixedVector<C2, SIZE>::const_iterator itv(v.begin());
    while (it != itend)
      *it++ += s * *itv++;
  }
   
  template <class C, unsigned int SIZE>
  template <class C2>
  void FixedVector<C, SIZE>::sadd(const C s, const FixedVector<C2, SIZE>& v)
  {
    assert(SIZE > 0);
    assert(SIZE == v.size());

    iterator it(begin()), itend(end());
    typename FixedVector<C2, SIZE>::const_iterator itv(v.begin());
    while(it != itend)
      {
        *it = s*(*it) + *itv++;
        ++it;
      }
  }
  
  template <class C, unsigned int SIZE>
  void FixedVector<C, SIZE>::scale(const C s)
  {
    assert(SIZE > 0);
    
    iterator it(begin()), itend(end());
    while(it != itend)
      *it++ *= s;
  }

  template <class C, unsigned int SIZE>
  template <class C2>
  inline
  FixedVector<C, SIZE>& FixedVector<C, SIZE>::operator += (const FixedVector<C2, SIZE>& v)
  {
    add(v); // handles the assertions
    return *this;
  }

  template <class C, unsigned int SIZE>
  template <class C2>
  void FixedVector<C, SIZE>::subtract(const FixedVector<C2, SIZE>& v)
  {
    assert(SIZE > 0);
    assert(SIZE == v.size());

    iterator it(begin()), itend(end());
    typename FixedVector<C2, SIZE>::const_iterator itv(v.begin());
    while (it != itend)
      *it++ -= *itv++;
  }

  template <class C, unsigned int SIZE>
  template <class C2>
  inline
  FixedVector<C, SIZE>& FixedVector<C, SIZE>::operator -= (const FixedVector<C2, SIZE>& v)
  {
    subtract(v); // handles the assertions
    return *this;
  }
   
  template <class C, unsigned int SIZE>
  FixedVector<C, SIZE>& FixedVector<C, SIZE>::operator *= (const C s)
  {
    scale(s); // handles the assertions
    return *this;
  }

  template <class C, unsigned int SIZE>
  FixedVector<C, SIZE>& FixedVector<C, SIZE>::operator /= (const C s)
  {
    // we don't catch the division by zero exception here!
    return (*this *= 1.0/s);
  }

  template <class C, unsigned int SIZE>
  inline
  std::ostream& operator << (std::ostream& os, const FixedVector<C, SIZE>& A)
  {
    print_vector(A, os);
    return os;
  }
  
}

// implementation of InfiniteVector inline functions

#include <cmath>
#include <algorithm>

namespace MathTL
{
  template <class C, class I>
  InfiniteVector<C,I>::InfiniteVector()
    : std::map<I,C>()
  {
  }

  template <class C, class I>
  InfiniteVector<C,I>::InfiniteVector(const InfiniteVector<C,I>& v)
    : std::map<I,C>(v)
  {
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::operator == (const InfiniteVector<C,I>& v) const
  {
    return std::equal(begin(), end(), v.begin());
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::operator != (const InfiniteVector<C,I>& v) const
  {
    return !((*this) == v);
  }

  template <class C, class I>
  inline
  C InfiniteVector<C,I>::operator [] (const I& index) const
  {
    // We must not use the map::operator [] for reading,
    // since it may add unwanted zero elements!
    typename std::map<I,C>::const_iterator it(lower_bound(index));
    if (it != std::map<I,C>::end() &&
	!std::map<I,C>::key_comp()(index,it->first))
      return it->second; 
    return C(0);
  }

  template <class C, class I>
  inline
  C& InfiniteVector<C,I>::operator [] (const I& index)
  {
    return std::map<I,C>::operator [] (index);
  }

  template <class C, class I>
  template <class C2>
  void InfiniteVector<C,I>::add(const InfiniteVector<C2,I>& v)
  {
    typename InfiniteVector<C2,I>::const_iterator itv(v.begin()), itvend(v.end());
    while (itv != itvend)
      this->operator [] (itv.index()) += *itv++;
  }
   
  template <class C, class I>
  template <class C2>
  void InfiniteVector<C,I>::add(const C2 s, const InfiniteVector<C2,I>& v)
  {
    // the following code can be optimized (not O(N) now)
    typename InfiniteVector<C2,I>::const_iterator itv(v.begin()), itvend(v.end());
    while (itv != itvend)
      this->operator [] (itv->index()) += s * *itv++;
  }
   
  template <class C, class I>
  template <class C2>
  void InfiniteVector<C,I>::sadd(const C s, const InfiniteVector<C2,I>& v)
  {
    // the following code can be optimized (not O(N) now)
    typename InfiniteVector<C2,I>::const_iterator itv(v.begin()), itvend(v.end());
    while (itv != itvend)
      this->operator [] (itv->index()) = 
	s*this->operator [] (itv->index()) + *itv++;
  }

  template <class C, class I>
  void InfiniteVector<C,I>::scale(const C s)
  {
    typename std::map<I,C>::iterator it(std::map<I,C>::begin()),
      itend(std::map<I,C>::end());
    while(it != itend)
      (*it++).second *= s;
  }

  template <class C, class I>
  template <class C2>
  inline
  InfiniteVector<C,I>& InfiniteVector<C,I>::operator += (const InfiniteVector<C2,I>& v)
  {
    add(v);
    return *this;
  }

  template <class C, class I>
  template <class C2>
  InfiniteVector<C,I>& InfiniteVector<C,I>::operator -= (const InfiniteVector<C2,I>& v)
  {
    typename InfiniteVector<C2,I>::const_iterator itv(v.begin()), itvend(v.end());
    while (itv != itvend)
      this->operator [] (itv.index()) -= *itv++;
    
    return *this;
  }
   
  template <class C, class I>
  InfiniteVector<C,I>& InfiniteVector<C,I>::operator *= (const C s)
  {
    scale(s);
    return *this;
  }

  template <class C, class I>
  InfiniteVector<C,I>& InfiniteVector<C,I>::operator /= (const C s)
  {
    // we don't catch the division by zero exception here!
    return (*this *= 1.0/s);
  }

//   template <class C, class I>
//   template <class C2>
//   const C InfiniteVector<C,I>::operator * (const InfiniteVector<C2,I>& v) const
//   {
//     if (this == reinterpret_cast<const InfiniteVector<C,I>*>(&v))
//       return l2_norm_sqr(*this);

//     C r(0);
    
//     const_iterator it(begin()), itend(end());
//     typename InfiniteVector<C2,I>::const_iterator itv(v.begin());
//     while(it != itend)
//       r += (*it++).value() * (*itv++).value();
    
//     return r;
//   }



  
//   template <class C, class I>
//   InfiniteVector<C,I>::Accessor::
//   Accessor(const typename std::map<I,C>::const_iterator& entry)
//     : entry_(entry)
//   {
//   }
  
//   template <class C, class I>
//   I InfiniteVector<C,I>::Accessor::index() const
//   {
//     return entry_->first;
//   }
  
//   template <class C, class I>
//   C InfiniteVector<C,I>::Accessor::value() const
//   {
//     return entry_->second;
//   }

//   template <class C, class I>
//   inline
//   bool
//   InfiniteVector<C,I>::Accessor::operator == (const Accessor& a) const
//   {
//     return (index() == a.index() && value() == a.value());
//   }

//   template <class C, class I>
//   inline
//   bool
//   InfiniteVector<C,I>::Accessor::operator != (const Accessor& a) const
//   {
//     return !((*this) == a);
//   }

  template <class C, class I>
  InfiniteVector<C,I>::const_iterator::
  const_iterator(const typename std::map<I,C>::const_iterator& entry)
    : std::map<I,C>::const_iterator(entry)
  {
  }

  template <class C, class I>
  typename InfiniteVector<C,I>::const_iterator
  InfiniteVector<C,I>::begin() const
  {
    return const_iterator(std::map<I,C>::begin());
  }

  template <class C, class I>
  typename InfiniteVector<C,I>::const_iterator
  InfiniteVector<C,I>::end() const
  {
    return const_iterator(std::map<I,C>::end());
  }

  template <class C, class I>
  inline
  const C&
  InfiniteVector<C,I>::const_iterator::operator * () const
  {
    return (std::map<I,C>::const_iterator::operator *()).second;
  }

  template <class C, class I>
  inline
  const C*
  InfiniteVector<C,I>::const_iterator::operator -> () const
  {
    return &(std::map<I,C>::const_iterator::operator *).second;
  }

  template <class C, class I>
  inline
  I InfiniteVector<C,I>::const_iterator::index() const
  {
    return (std::map<I,C>::const_iterator::operator *()).first;
  }

  template <class C, class I>
  inline
  typename InfiniteVector<C,I>::const_iterator&
  InfiniteVector<C,I>::const_iterator::operator ++ ()
  {
    std::map<I,C>::const_iterator::operator ++ ();
    return *this;
  }

  template <class C, class I>
  inline
  typename InfiniteVector<C,I>::const_iterator
  InfiniteVector<C,I>::const_iterator::operator ++ (int step)
  {
    InfiniteVector<C,I>::const_iterator r(*this);
    std::map<I,C>::const_iterator::operator ++ (step);
    return r;
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::const_iterator::
  operator == (const const_iterator& it) const
  {
    return (index() == it.index());
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::const_iterator::
  operator != (const const_iterator& it) const
  {
    return !(*this == it);
  }

  template <class C, class I>
  std::ostream& operator << (std::ostream& os,
			     const InfiniteVector<C,I>& v)
  {
    if (v.begin() ==  v.end())
      {
	os << "0";
      }
    else
      {
	for (typename InfiniteVector<C,I>::const_iterator it(v.begin());
	     it != v.end(); ++it)
	  {
	    os << it.index() << ": " << *it << endl;
	  }
      }

    return os;
  }
}

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
  InfiniteVector<C,I>::Accessor::
  Accessor(const typename std::map<I,C>::const_iterator& entry)
    : entry_(entry)
  {
  }
  
  template <class C, class I>
  I InfiniteVector<C,I>::Accessor::index() const
  {
    return entry_->first;
  }
  
  template <class C, class I>
  C InfiniteVector<C,I>::Accessor::value() const
  {
    return entry_->second;
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::Accessor::operator == (const Accessor& a) const
  {
    return (index() == a.index() && value() == a.value());
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::Accessor::operator != (const Accessor& a) const
  {
    return !((*this) == a);
  }

  template <class C, class I>
  InfiniteVector<C,I>::const_iterator::
  const_iterator(const typename std::map<I,C>::const_iterator& entry)
    : accessor_(entry)
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
  const typename InfiniteVector<C,I>::Accessor&
  InfiniteVector<C,I>::const_iterator::operator * () const
  {
    return accessor_;
  }

  template <class C, class I>
  inline
  const typename InfiniteVector<C,I>::Accessor*
  InfiniteVector<C,I>::const_iterator::operator -> () const
  {
    return &accessor_;
  }

  template <class C, class I>
  inline
  typename InfiniteVector<C,I>::const_iterator&
  InfiniteVector<C,I>::const_iterator::operator ++ ()
  {
    ++accessor_.entry_;
    return *this;
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::const_iterator::
  operator == (const const_iterator& it) const
  {
    return (accessor_.index() == it.accessor_.index());
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
	    os << it->index() << ": " << it->value() << endl;
	  }
      }

    return os;
  }
}

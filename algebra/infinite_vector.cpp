// implementation of InfiniteVector inline functions

#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

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
  InfiniteVector<C,I>&
  InfiniteVector<C,I>::operator = (const InfiniteVector<C,I>& v)
  {
    std::map<I,C>::operator = (v);
    return *this;
  }

  template <class C, class I>
  inline
  void InfiniteVector<C,I>::swap(InfiniteVector<C,I>& v)
  {
    std::map<I,C>::swap(v);
  }

  template <class C, class I>
  inline
  void InfiniteVector<C,I>::clear()
  {
    std::map<I,C>::clear();
  }

  template <class C, class I>
  inline
  size_t InfiniteVector<C,I>::size() const
  {
    return std::map<I,C>::size();
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

  template <class C, class I>
  template <class C2>
  const C InfiniteVector<C,I>::operator * (const InfiniteVector<C2,I>& v) const
  {
    if (this == reinterpret_cast<const InfiniteVector<C,I>*>(&v))
      return l2_norm_sqr(*this);

    C r(0);
    
    typename InfiniteVector<C,I>::const_iterator it(begin()), itend(end()),
      itv(v.begin()), itvend(v.end());
    for (; it != itend && itv != itvend; ++it)
      {
 	while (itv != itvend && itv < it) ++itv;
 	if (it == itv)
 	  r += *it * *itv;
      }

    return r;
  }

  template <class C, class I>
  double operator * (const InfiniteVector<C,I>& v, const InfiniteVector<C,I>& w)
  {
    double r(0);
    typedef typename InfiniteVector<C,I>::entry_const_iterator entry_const_iterator;
    entry_const_iterator itv(v.begin()), itvend(v.end()), itw(w.begin()), itwend(w.end());
    for (; itv != itvend && itw != itwend; ++itv)
      {
 	while (itw != itwend && itw.index() < itv.index()) ++itw;
 	if (itv.index() == itw.index())
 	  r += itw.entry() * itv.entry();
      }
    return r;
  }

  template <class C, class I>
  double InfiniteVector<C,I>::weak_norm(const double tau) const
  {
    double r(0.0);

    if (size() > 0)
      {
	// prepare vector to be sorted
	std::vector<std::pair<I,C> > sv(size());
	unsigned int id(0);
	for (const_iterator it(begin()), itend(end());
	     it != itend; ++it, ++id)
	  {
	    sv[id] = std::make_pair<I,C>(it.index(), *it);
	  }
	  
	// sort vector (Introsort, O(N*log N))
	sort(sv.begin(), sv.end(), decreasing_order());
	  
	// compute \|*this\|_{\ell^w_\tau}:=\sup_{N=1}^\infty N^{1/tau}|v_N^*|
	// where the v_N^* are the decreasing rearrangement of me
	for (unsigned int N(1); N <= sv.size(); N++)
	  r = std::max(r, pow(N, 1.0/tau) * fabs(sv[N-1].second));
      }

    return r;
  }

  template <class C, class I>
  void InfiniteVector<C,I>::n_coarse(const double eps, InfiniteVector<C,I>& v) const
  {
    // We use a straightforward implementation with complexity O(N*log(N)):
    // - sort my entries in modulus
    //   1. possibility: use a helper multimap object
    //   2. possibility: use a sorted vector (preferred solution, since no
    //                   slow insertion sort algorithm is launched!)
    // - insert the largest in modulus entries into v until
    //     \|*this-v\|_{\ell_2}\le\epsilon
    //
    // Another possibility would be binary binning, which we will implement
    // in a later stage of the library!

    v.clear();
    if (size() > 0)
      {
	// prepare vector to be sorted
	std::vector<std::pair<I,C> > sv(size());
	unsigned int id(0);
	for (const_iterator it(begin()), itend(end());
	     it != itend; ++it, ++id)
	  {
	    sv[id] = std::make_pair<I,C>(it.index(), *it);
	  }
	  
	// sort vector (Introsort, O(N*log N))
	sort(sv.begin(), sv.end(), decreasing_order());

	// insert largest in modulus entries until tolerance is reached
	double coarsenorm(0);
	double nrm(l2_norm(*this));
	double bound(nrm*nrm - eps*eps);
	typename std::vector<std::pair<I,C> >::iterator it(sv.begin());
	do
	  {
	    coarsenorm += it->second * it->second;
	    ++it;
	  }
	while ((it != sv.end()) && (coarsenorm < bound));
	sv.erase(it, sv.end());

	// insert relevant entries in v (-> insertion sort, we hope that
	// the number of entries is neglectible)
	for (unsigned int i(0); i < sv.size(); i++)
	  v[sv[i].first] = sv[i].second;
      }
  }

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
    // quick, dirty hack
    return (static_cast<typename std::map<I,C>::const_iterator>(*this)
	    == static_cast<typename std::map<I,C>::const_iterator>(it));
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
  inline
  bool
  InfiniteVector<C,I>::const_iterator::
  operator < (const const_iterator& it) const
  {
    return (index() < it.index());
  }

  template <class C, class I>
  inline
  void swap(InfiniteVector<C,I>& v1, InfiniteVector<C,I>& v2)
  {
    v1.swap(v2);
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

// implementation of InfiniteVector inline functions

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>
#include <utils/tiny_tools.h>

using std::cout;
using std::endl;

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
    return get_coefficient(index);
  }

  template <class C, class I>
  C InfiniteVector<C,I>::get_coefficient(const I& index) const
  {
    typename std::map<I,C>::const_iterator it(this->lower_bound(index));
    if (it != std::map<I,C>::end())
      {
	if (!std::map<I,C>::key_comp()(index, it->first))
	  return it->second;
      }
    return C(0);
  }

  template <class C, class I>
  C& InfiniteVector<C,I>::operator [] (const I& index)
  {
    // efficient add-or-update, cf. Meyers, Effective STL
    typename std::map<I,C>::iterator it(this->lower_bound(index));
    if (it != std::map<I,C>::end() &&
	!std::map<I,C>::key_comp()(index, it->first))
      return it->second;
    else
      return std::map<I,C>::insert(it, typename std::map<I,C>::value_type(index, C(0)))->second;
// alternative: return std::map<I,C>::insert(typename std::map<I,C>::value_type(index, C(0))).first->second;
  }

  template <class C, class I>
  void InfiniteVector<C,I>::set_coefficient(const I& index, const C value)
  {
    typename std::map<I,C>::iterator it(this->lower_bound(index));
    if (it != std::map<I,C>::end() &&
	!std::map<I,C>::key_comp()(index,it->first)) {
      it->second = value;
    }
    else {
      std::map<I,C>::insert(it, typename std::map<I,C>::value_type(index, value));
      //std::map<I,C>::insert(typename std::map<I,C>::value_type(index, value));
    }
  }

  template <class C, class I>
  void InfiniteVector<C,I>::add_coefficient(const I& index, const C increment)
  {
    // efficient add-or-update, cf. Meyers, Effective STL
    typename std::map<I,C>::iterator it(this->lower_bound(index));
    if (it != std::map<I,C>::end() &&
	!std::map<I,C>::key_comp()(index, it->first)) {
      // we already have a nontrivial coefficient
      if ((it->second += increment) == C(0))
	std::map<I,C>::erase(it);
    } else {
      // insert the increment as new coefficient
      std::map<I,C>::insert(it, typename std::map<I,C>::value_type(index, increment));
    }
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
  void
  InfiniteVector<C,I>::support(std::set<I>& supp) const
  {
    supp.clear();
    for (const_iterator it(begin()), itend(end()); it != itend; ++it)
      supp.insert(supp.end(), it.index());
  }

  template <class C, class I>
  void
  InfiniteVector<C,I>::clip(const std::set<I>& supp)
  {
    std::map<I,C> v;

    const_iterator it(begin()), itend(end());
    typename std::set<I>::const_iterator suppit(supp.begin()), suppend(supp.end());
    for (; it != itend && suppit != suppend; ++it)
      {
	while (suppit != suppend && *suppit < it.index()) ++suppit;
	if (*suppit == it.index())
	  v.insert(v.end(), std::pair<I,C>(it.index(), *it));
      }

    std::map<I,C>::swap(v);
  }

  template <class C, class I>
  void InfiniteVector<C,I>::compress(const double eta)
  {
    // a hardcore STL implementation inspired by Meyers, Effective STL:
    std::map<I,C> v;
    remove_copy_if(std::map<I,C>::begin(),
		   std::map<I,C>::end(),
		   std::inserter(v, v.end()),
		   threshold_criterion<I,C>(eta));
    std::map<I,C>::swap(v);
  }

  template <class C, class I>
  void InfiniteVector<C,I>::add(const InfiniteVector<C,I>& v)
  {
#if 1
    std::map<I,C> help;

    // The following O(N) algorithm is adapted from the STL algorithm set_union(),
    // cf. stl_algo.h ...

    typename InfiniteVector<C,I>::const_iterator it(begin()), itend(end()),
      itv(v.begin()), itvend(v.end());
    typename std::map<I,C>::iterator hint(help.begin()), hint2(help.begin());

    while (it != itend && itv != itvend)
      {
	if (it.index() < itv.index())
	  {
	    hint2 = help.insert(hint, std::pair<I,C>(it.index(), *it));
	    hint = hint2;
	    ++it;
	  }
	else
	  {
	    if (itv.index() < it.index())
	      {
		hint2 = help.insert(hint, std::pair<I,C>(itv.index(), *itv));
		hint = hint2;
		++itv;
	      }
	    else
	      {
		const C value(*it + *itv);
		if (value != C(0)) {
		  hint2 = help.insert(hint, std::pair<I,C>(itv.index(), value));
		  hint = hint2;
		}
		++it;
		++itv;
	      }
	  }
      }

    while (it != itend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(it.index(), *it));
	hint = hint2;
	++it;
      }

    while (itv != itvend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(itv.index(), *itv));
	hint = hint2;
	++itv;
      }

    std::map<I,C>::swap(help);
#else
    // the following code can be optimized (not O(N) now)
    typename InfiniteVector<C,I>::const_iterator itv(v.begin()), itvend(v.end());
    for (; itv != itvend; ++itv)
      {
	C help(get_coefficient(itv.index()) + *itv);

	if (help != C(0))
	  set_coefficient(itv.index(), help);
	else
	  std::map<I,C>::erase(itv.index());
      }
#endif
  }

  template <class C, class I>
  void InfiniteVector<C,I>::add(const C s, const InfiniteVector<C,I>& v)
  {
#if 1
    std::map<I,C> help;

    // The following O(N) algorithm is adapted from the STL algorithm set_union(),
    // cf. stl_algo.h ...

    typename InfiniteVector<C,I>::const_iterator it(begin()), itend(end()),
      itv(v.begin()), itvend(v.end());
    typename std::map<I,C>::iterator hint(help.begin()), hint2(help.begin());

    while (it != itend && itv != itvend)
      {
	if (it.index() < itv.index())
	  {
	    hint2 = help.insert(hint, std::pair<I,C>(it.index(), *it));
	    hint = hint2;
	    ++it;
	  }
	else
	  {
	    if (itv.index() < it.index())
	      {
		hint2 = help.insert(hint, std::pair<I,C>(itv.index(), s * *itv));
		hint = hint2;
		++itv;
	      }
	    else
	      {
		const C value(*it + s * *itv);
		if (value != C(0)) {
		  hint2 = help.insert(hint, std::pair<I,C>(itv.index(), value));
		  hint = hint2;
		}
		++it;
		++itv;
	      }
	  }
      }

    while (it != itend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(it.index(), *it));
	hint = hint2;
	++it;
      }

    while (itv != itvend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(itv.index(), s * *itv));
	hint = hint2;
	++itv;
      }

    std::map<I,C>::swap(help);
#else
    // the following code can be optimized (not O(N) now)
    typename InfiniteVector<C,I>::const_iterator itv(v.begin()), itvend(v.end());
    for (; itv != itvend; ++itv)
      {
	C help(get_coefficient(itv.index()) + s * *itv);
	if (help != C(0))
	  set_coefficient(itv.index(), help);
	else
	  std::map<I,C>::erase(itv.index());
      }
#endif
  }

  template <class C, class I>
  void InfiniteVector<C,I>::sadd(const C s, const InfiniteVector<C,I>& v)
  {
#if 1
    std::map<I,C> help;

    // The following O(N) algorithm is adapted from the STL algorithm set_union(),
    // cf. stl_algo.h ...

    typename InfiniteVector<C,I>::const_iterator it(begin()), itend(end()),
      itv(v.begin()), itvend(v.end());
    typename std::map<I,C>::iterator hint(help.begin()), hint2(help.begin());

    while (it != itend && itv != itvend)
      {
	if (it.index() < itv.index())
	  {
	    hint2 = help.insert(hint, std::pair<I,C>(it.index(), s * *it));
	    hint = hint2;
	    ++it;
	  }
	else
	  {
	    if (itv.index() < it.index())
	      {
		hint2 = help.insert(hint, std::pair<I,C>(itv.index(), *itv));
		hint = hint2;
		++itv;
	      }
	    else
	      {
		const C value(s * *it + *itv);
		if (value != C(0)) {
		  hint2 = help.insert(hint, std::pair<I,C>(itv.index(), value));
		  hint = hint2;
		}
		++it;
		++itv;
	      }
	  }
      }

    while (it != itend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(it.index(), s * *it));
	hint = hint2;
	++it;
      }

    while (itv != itvend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(itv.index(), *itv));
	hint = hint2;
	++itv;
      }

    std::map<I,C>::swap(help);
#else
    // the following code can be optimized (not O(N) now)
    typename InfiniteVector<C,I>::const_iterator itv(v.begin()), itvend(v.end());
    for (; itv != itvend; ++itv)
      {
	C help(s * get_coefficient(itv.index()) + *itv);
	if (help != C(0))
	  set_coefficient(itv.index(), help);
	else
	  std::map<I,C>::erase(itv.index());
      }
#endif
  }

  template <class C, class I>
  void InfiniteVector<C,I>::scale(const C s)
  {
    if (s == C(0))
      clear();
    else
      {
	typename std::map<I,C>::iterator it(std::map<I,C>::begin()),
	  itend(std::map<I,C>::end());
	while(it != itend)
	  (*it++).second *= s;
      }
  }

  template <class C, class I>
  void InfiniteVector<C,I>::scale(const InfiniteDiagonalMatrix<C,I>* D, const int k)
  {
    for (typename std::map<I,C>::iterator it(std::map<I,C>::begin()),
	   itend(std::map<I,C>::end());
	 it != itend; ++it)
      it->second *= pow(D->diag(it->first), k);
  }

  template <class C, class I>
  inline
  InfiniteVector<C,I>& InfiniteVector<C,I>::operator += (const InfiniteVector<C,I>& v)
  {
    add(v);
    return *this;
  }

  template <class C, class I>
  void InfiniteVector<C,I>::subtract(const InfiniteVector<C,I>& v)
  {
#if 1
    std::map<I,C> help;

    // The following O(N) algorithm is adapted from the STL algorithm set_union(),
    // cf. stl_algo.h ...

    typename InfiniteVector<C,I>::const_iterator it(begin()), itend(end()),
      itv(v.begin()), itvend(v.end());
    typename std::map<I,C>::iterator hint(help.begin()), hint2(help.begin());

    while (it != itend && itv != itvend)
      {
	if (it.index() < itv.index())
	  {
	    hint2 = help.insert(hint, std::pair<I,C>(it.index(), *it));
	    hint = hint2;
	    ++it;
	  }
	else
	  {
	    if (itv.index() < it.index())
	      {
		hint2 = help.insert(hint, std::pair<I,C>(itv.index(), - *itv));
		hint = hint2;
		++itv;
	      }
	    else
	      {
		const C value(*it - *itv);
		if (value != C(0)) {
		  hint2 = help.insert(hint, std::pair<I,C>(itv.index(), value));
		  hint = hint2;
		}
		++it;
		++itv;
	      }
	  }
      }

    while (it != itend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(it.index(), *it));
	hint = hint2;
	++it;
      }

    while (itv != itvend)
      {
	hint2 = help.insert(hint, std::pair<I,C>(itv.index(), - *itv));
	hint = hint2;
	++itv;
      }

    std::map<I,C>::swap(help);
#else
    // the following code can be optimized (not O(N) now)
    typename InfiniteVector<C,I>::const_iterator itv(v.begin()), itvend(v.end());
    for (; itv != itvend; ++itv)
      {
	C help(get_coefficient(itv.index()) - *itv);
	if (help != C(0))
	  set_coefficient(itv.index(), help);
	else
	  std::map<I,C>::erase(itv.index());
      }
#endif
  }

  template <class C, class I>
  inline
  InfiniteVector<C,I>& InfiniteVector<C,I>::operator -= (const InfiniteVector<C,I>& v)
  {
    subtract(v);
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
  const C InfiniteVector<C,I>::operator * (const InfiniteVector<C,I>& v) const
  {
    if (this == reinterpret_cast<const InfiniteVector<C,I>*>(&v))
      return l2_norm_sqr(*this);

    C r(0);
    
    for (typename InfiniteVector<C,I>::const_iterator
	   it(begin()),
	   itend(end()),
	   itv(v.begin()),
	   itvend(v.end());
	 it != itend && itv != itvend; ++it)
      {
 	while (itv != itvend && itv.index() < it.index()) ++itv;
	if (itv != itvend)
	  if (it.index() == itv.index())
	    r += *it * *itv;
      }

    return r;
  }

  template <class C, class I>
  double operator * (const InfiniteVector<C,I>& v, const InfiniteVector<C,I>& w)
  {
    double r(0);
    typedef typename InfiniteVector<C,I>::const_iterator const_iterator;
    const_iterator itv(v.begin()), itvend(v.end()), itw(w.begin()), itwend(w.end());
    for (; itv != itvend && itw != itwend; ++itv)
      {
 	while (itw != itwend && itw.index() < itv.index()) ++itw;
 	if (itv.index() == itw.index())
 	  r += *itw * *itv;
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
	    sv[id] = std::pair<I, C>(it.index(), *it); // can't use make_pair for gcc 2.95
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
  void InfiniteVector<C,I>::COARSE(const double eps, InfiniteVector<C,I>& v) const
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
    if (size() > 0) {
      if (eps == 0)
	v = *this;
      else {
	
	// prepare vector to be sorted
	
	std::vector<std::pair<I,C> > sv(size());
	

	//abort();
	unsigned int id(0);
	for (const_iterator it(begin()), itend(end()); it != itend; ++it, ++id) {
	  //cout << it.index() << endl;
	  sv[id] = std::pair<I,C>(it.index(), *it); // can't use make_pair for gcc 2.95
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
	  v.set_coefficient(sv[i].first, sv[i].second);
      }
    }
  }

  template <class C, class I>
  const double
  InfiniteVector<C,I>::wrmsqr_norm(const double atol, const double rtol,
				   const InfiniteVector<C,I>& v, const InfiniteVector<C,I>& w) const
  {
    double result = 0;
    
    for (const_iterator it(begin()), itend(end());
 	 it != itend; ++it)
      {
  	const double help = *it / (atol + rtol * std::max(fabs(v.get_coefficient(it.index())),
							  fabs(w.get_coefficient(it.index()))));
  	result += help * help;
      }

    return result == 0 ? 0 : sqrt(result/size());
  }

  template <class C, class I>
  inline
  InfiniteVector<C,I>::const_iterator::
  const_iterator(const typename std::map<I,C>::const_iterator& entry)
    : std::map<I,C>::const_iterator(entry)
  {
  }

  template <class C, class I>
  inline
  typename InfiniteVector<C,I>::const_iterator
  InfiniteVector<C,I>::begin() const
  {
    return const_iterator(std::map<I,C>::begin());
  }

  template <class C, class I>
  inline
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
    return &((std::map<I,C>::const_iterator::operator *()).second);
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
    typename InfiniteVector<C,I>::const_iterator r(*this);
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
  InfiniteVector<C,I>::const_reverse_iterator::
  const_reverse_iterator(const std::reverse_iterator<typename std::map<I,C>::const_iterator>& entry)
    : std::reverse_iterator<typename std::map<I,C>::const_iterator>(entry)
  {
  }

  template <class C, class I>
  typename InfiniteVector<C,I>::const_reverse_iterator
  InfiniteVector<C,I>::rbegin() const
  {
    return const_reverse_iterator(std::reverse_iterator<typename std::map<I,C>::const_iterator>
				  (std::map<I,C>::end()));
  }

  template <class C, class I>
  typename InfiniteVector<C,I>::const_reverse_iterator
  InfiniteVector<C,I>::rend() const
  {
    return const_reverse_iterator(std::reverse_iterator<typename std::map<I,C>::const_iterator>
				  (std::map<I,C>::begin()));
  }

  template <class C, class I>
  inline
  const C&
  InfiniteVector<C,I>::const_reverse_iterator::operator * () const
  {
    return (std::reverse_iterator<typename std::map<I,C>::const_iterator>::operator *()).second;
  }

  template <class C, class I>
  inline
  const C*
  InfiniteVector<C,I>::const_reverse_iterator::operator -> () const
  {
    return &(std::reverse_iterator<typename std::map<I,C>::const_iterator>::operator *).second;
  }

  template <class C, class I>
  inline
  I InfiniteVector<C,I>::const_reverse_iterator::index() const
  {
    return (std::reverse_iterator<typename std::map<I,C>::const_iterator>::operator *()).first;
  }

  template <class C, class I>
  inline
  typename InfiniteVector<C,I>::const_reverse_iterator&
  InfiniteVector<C,I>::const_reverse_iterator::operator ++ ()
  {
    std::reverse_iterator<typename std::map<I,C>::const_iterator>::operator ++ ();
    return *this;
  }

  template <class C, class I>
  inline
  typename InfiniteVector<C,I>::const_reverse_iterator
  InfiniteVector<C,I>::const_reverse_iterator::operator ++ (int step)
  {
    typename InfiniteVector<C,I>::const_reverse_iterator r(*this);
    std::reverse_iterator<typename std::map<I,C>::const_iterator>::operator ++ (step);
    return r;
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::const_reverse_iterator::
  operator == (const const_reverse_iterator& it) const
  {
    // quick, dirty hack
    return (static_cast<std::reverse_iterator<typename std::map<I,C>::const_iterator> >(*this)
	    == static_cast<std::reverse_iterator<typename std::map<I,C>::const_iterator> >(it));
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::const_reverse_iterator::
  operator != (const const_reverse_iterator& it) const
  {
    return !(*this == it);
  }

  template <class C, class I>
  inline
  bool
  InfiniteVector<C,I>::const_reverse_iterator::
  operator < (const const_reverse_iterator& it) const
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
	    os << it.index() << ": " << *it << std::endl;
	  }
      }

    return os;
  }
}

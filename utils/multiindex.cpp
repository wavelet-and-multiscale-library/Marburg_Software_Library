// implementation for multiindex.h
#include <cassert>
#include <algorithm>

namespace MathTL
{
  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex()
    : FixedArray1D<I, DIMENSION>()
  {
    for (unsigned int i(0); i < DIMENSION; i++)
      FixedArray1D<I, DIMENSION>::operator [] (i) = I(0);
  }

  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex(const MultiIndex<I, DIMENSION>& lambda)
    : FixedArray1D<I, DIMENSION>(lambda)
  {
  }

  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex(const I& i0)
    : FixedArray1D<I, DIMENSION>()
  {
    assert(DIMENSION == 1);
    FixedArray1D<I, DIMENSION>::operator [] (0) = i0;
  }

  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex(const I& i0, const I& i1)
    : FixedArray1D<I, DIMENSION>()
  {
    assert(DIMENSION == 2);
    FixedArray1D<I, DIMENSION>::operator [] (0) = i0;
    FixedArray1D<I, DIMENSION>::operator [] (1) = i1;
  }

  template <class I, unsigned int DIMENSION>
  bool MultiIndex<I, DIMENSION>::operator == (const MultiIndex& lambda) const
  {
    for (unsigned int i(0); i < DIMENSION; i++)
      if (FixedArray1D<I, DIMENSION>::operator [] (i) != lambda[i]) return false;
    return true;
  }

  template <class I, unsigned int DIMENSION>
  inline
  bool MultiIndex<I, DIMENSION>::operator < (const MultiIndex& lambda) const
  {
    return std::lexicographical_compare(FixedArray1D<I, DIMENSION>::begin(),
					FixedArray1D<I, DIMENSION>::end(),
					lambda.begin(), lambda.end());
  }

  template<class I, unsigned int DIMENSION>
  std::set<MultiIndex<I, DIMENSION> >
  cuboid_indices(const MultiIndex<I, DIMENSION>& alpha,
		 const MultiIndex<I, DIMENSION>& beta)
  {
    typedef std::set<MultiIndex<I, DIMENSION> > set_type;
    set_type r;

    MultiIndex<I, DIMENSION> gamma(alpha);
    for (I index(alpha[0]); index <= beta[0]; index++)
      {
	gamma[0] = index;
	r.insert(gamma);
      }
    for (unsigned int i(1); i < DIMENSION; i++)
      {
	set_type sofar(r);
	for (typename set_type::const_iterator it(sofar.begin()); it != sofar.end(); ++it)
	  {
	    gamma = *it;
	    for (I index(alpha[i]); index <= beta[i]; index++)
	      {
		gamma[i] = index;
		r.insert(gamma);
	      }
	  }
      }

    return r;
  }
}

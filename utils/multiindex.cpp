// implementation for multiindex.h
#include <cassert>
#include <algorithm>
#include <utils/tiny_tools.h>

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

  template <unsigned int DIMENSION>
  std::set<MultiIndex<unsigned int, DIMENSION> >
  degree_indices(const unsigned int k)
  {
    typedef std::set<MultiIndex<unsigned int, DIMENSION> > set_type;
    set_type r;

    MultiIndex<unsigned int, DIMENSION> alpha;
    r.insert(alpha);
    for (unsigned int i(1); i <= k; i++)
      {
	set_type sofar(r);
	r.clear();
	for (typename set_type::const_iterator it(sofar.begin()); it != sofar.end(); ++it)
	  {
	    alpha = *it;
	    for (unsigned int j(0); j < DIMENSION; j++)
	      {
		alpha[j] += 1;
		r.insert(alpha);
		alpha[j] -= 1;
	      }
	  }	
      }

    return r;
  }

  template <unsigned int DIMENSION>
  unsigned int multi_degree(const MultiIndex<unsigned int, DIMENSION>& alpha)
  {
    unsigned int r(alpha[0]);
    for (unsigned int i(1); i < DIMENSION; i++)
      r += alpha[i];
    return r;
  }

  template<unsigned int DIMENSION>
  unsigned int multi_faculty(const MultiIndex<unsigned int, DIMENSION>& alpha)
  {
    unsigned int r(faculty(alpha[0]));
    for (unsigned int i(1); i < DIMENSION; i++)
      r *= faculty(alpha[i]);
    return r;
  }

  template <class I, unsigned int DIMENSION>
  int multi_power(const MultiIndex<I, DIMENSION>& alpha,
		  const MultiIndex<unsigned int, DIMENSION>& beta)
  {
    int r(intpower(alpha[0], beta[0]));
    for (unsigned int i(1); i < DIMENSION; i++)
      r *= intpower(alpha[i], beta[i]);
    return r;
  }

  template <class I, unsigned int DIMENSION>
  int multi_binomial(const MultiIndex<I, DIMENSION>& beta,
		     const MultiIndex<I, DIMENSION>& alpha)
  {
    int r(binomial(beta[0], alpha[0]));
    for (unsigned int i(1); i < DIMENSION; i++)
      r *= binomial(beta[i], alpha[i]);
    return r;
  }
}

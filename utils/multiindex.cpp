// implementation for multiindex.h
#include <algorithm>

namespace MathTL
{
  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex()
    : FixedArray1D<I, DIMENSION>()
  {
  }

  template <class I, unsigned int DIMENSION>
  MultiIndex<I, DIMENSION>::MultiIndex(const MultiIndex<I, DIMENSION>& lambda)
    : FixedArray1D<I, DIMENSION>(lambda)
  {
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
    return std::lexicographical_compare(begin(), end(), lambda.begin(), lambda.end());
  }
}

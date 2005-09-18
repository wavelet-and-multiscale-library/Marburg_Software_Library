// implementation for tp_index.h

#include <cassert>

namespace WaveletTL
{
  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>::TensorProductIndex(const TensorProductBasis<BASIS1,BASIS2>* basis)
    : basis_(basis),
      index1_(basis == 0 ? 0 : &basis->basis1()),
      index2_(basis == 0 ? 0 : &basis->basis2())
  {
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>::TensorProductIndex(const TensorProductIndex<BASIS1,BASIS2>& lambda)
    : basis_(lambda.basis_), index1_(lambda.index1_), index2_(lambda.index2_)
  {
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>::TensorProductIndex(const TensorProductBasis<BASIS1,BASIS2>* basis,
							const typename BASIS1::Index& index1,
							const typename BASIS2::Index& index2)
    : basis_(basis), index1_(index1), index2_(index2)
  {
  }
  
  template <class BASIS1, class BASIS2>
  bool
  TensorProductIndex<BASIS1,BASIS2>::operator == (const TensorProductIndex& lambda) const
  {
    return index1_ == lambda.index1() && index2_ == lambda.index2();
  }

  template <class BASIS1, class BASIS2>
  bool
  TensorProductIndex<BASIS1,BASIS2>::operator < (const TensorProductIndex& lambda) const
  {
    // We want to have (j,e,k) < (j',e',k') iff
    //   j<j' or (j=j' and (e<e' or (e=e' and k<k')),
    // where e and k are lexicographically ordered vectors
    return (j() < lambda.j() ||
	    (j() == lambda.j() && ((index1().e() < lambda.index1().e() ||
				    (index1().e() == lambda.index1().e() && index2().e() < lambda.index2().e())) ||
				   ((index1().e() == lambda.index1().e() && index2().e() == lambda.index2().e()) &&
				    (index1().k() < lambda.index1().k() ||
				     (index1().k() == lambda.index1().k() && index2().k() < lambda.index2().k()))))));
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  first_generator(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j)
  {
    assert(j >= basis->j0());
    return TensorProductIndex<BASIS1,BASIS2>(basis,
					     first_generator(&basis->basis1(), j),
					     first_generator(&basis->basis2(), j));
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  last_generator(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j)
  {
    assert(j >= basis->j0());
    return TensorProductIndex<BASIS1,BASIS2>(basis,
					     last_generator(&basis->basis1(), j),
					     last_generator(&basis->basis2(), j));
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  first_wavelet(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j)
  {
    assert(j >= basis->j0());
    return TensorProductIndex<BASIS1,BASIS2>(basis,
					     first_generator(&basis->basis1(), j),
					     first_wavelet(&basis->basis2(), j));
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  last_wavelet(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j)
  {
    assert(j >= basis->j0());
    return TensorProductIndex<BASIS1,BASIS2>(basis,
					     last_wavelet(&basis->basis1(), j),
					     last_wavelet(&basis->basis2(), j));
  }
}

// implementation for tp_index.h

#include <cassert>

namespace WaveletTL
{
  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>::TensorProductIndex()
  {
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>::TensorProductIndex(const TensorProductIndex<BASIS1,BASIS2>& lambda)
    : index1_(lambda.index1_), index2_(lambda.index2_)
  {
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>::TensorProductIndex(const typename BASIS1::Index& index1,
							const typename BASIS2::Index& index2)
    : index1_(index1), index2_(index2)
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
    // where e and k are already lexicographically ordered vectors, respectively
    return (j() < lambda.j() ||
	    (j() == lambda.j() && ((index1().e() < lambda.index1().e() ||
				    (index1().e() == lambda.index1().e() && index2().e() < lambda.index2().e())) ||
				   ((index1().e() == lambda.index1().e() && index2().e() == lambda.index2().e()) &&
				    (index1().k() < lambda.index1().k() ||
				     (index1().k() == lambda.index1().k() && index2().k() < lambda.index2().k()))))));
  }

  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>&
  TensorProductIndex<BASIS1,BASIS2>::operator ++ ()
  {
    if (index2() == BASIS2::last_index(j(), index2().e())) {
      if (index1() == BASIS1::last_index(j(), index1().e())) {
	if (index2() == BASIS2::last_wavelet(j())) {
 	  if (index1() == BASIS1::last_wavelet(j())) {
 	    index1_ = BASIS1::first_generator(j()+1); // increments j
	    index2_ = BASIS2::first_wavelet(j()); // no generators on higher scales
	  } else {
	    ++index1_;
	    index2_ = BASIS2::first_generator(j());
	  }
	} else {
	  ++index2_;
	  index1_ = BASIS1::first_index(j(), index1().e());
	}
      } else {
	++index1_;
	index2_ = BASIS2::first_index(j(), index2().e());
      }
    } else
      ++index2_;
    
    return *this;
  }

//   template <class BASIS1, class BASIS2>
//   TensorProductIndex<BASIS1,BASIS2>
//   first_generator(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS1,BASIS2>(first_generator(&basis->basis1(), j),
// 					     first_generator(&basis->basis2(), j));
//   }

//   template <class BASIS1, class BASIS2>
//   TensorProductIndex<BASIS1,BASIS2>
//   last_generator(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS1,BASIS2>(last_generator(&basis->basis1(), j),
// 					     last_generator(&basis->basis2(), j));
//   }

//   template <class BASIS1, class BASIS2>
//   TensorProductIndex<BASIS1,BASIS2>
//   first_wavelet(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS1,BASIS2>(basis,
// 					     first_generator(&basis->basis1(), j),
// 					     first_wavelet(&basis->basis2(), j));
//   }

//   template <class BASIS1, class BASIS2>
//   TensorProductIndex<BASIS1,BASIS2>
//   last_wavelet(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS1,BASIS2>(last_wavelet(&basis->basis1(), j),
// 					     last_wavelet(&basis->basis2(), j));
//   }
}

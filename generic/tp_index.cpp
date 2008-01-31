// implementation for tp_index.h

#include <cassert>

namespace WaveletTL
{
  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>::TensorProductIndex()
  {
  }

  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>::TensorProductIndex(const TensorProductIndex<BASIS0,BASIS1>& lambda)
    : index0_(lambda.index0_), index1_(lambda.index1_)
  {
  }

  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>::TensorProductIndex(const typename BASIS0::Index& index0,
							const typename BASIS1::Index& index1)
    : index0_(index0), index1_(index1)
  {
  }
  
  template <class BASIS0, class BASIS1>
  bool
  TensorProductIndex<BASIS0,BASIS1>::operator == (const TensorProductIndex& lambda) const
  {
    return index0_ == lambda.index0() && index1_ == lambda.index1();
  }

  template <class BASIS0, class BASIS1>
  bool
  TensorProductIndex<BASIS0,BASIS1>::operator < (const TensorProductIndex& lambda) const
  {
    // We want to have (j,e,k) < (j',e',k') iff
    //   j<j' or (j=j' and (e<e' or (e=e' and k<k')),
    // where e and k are already lexicographically ordered vectors, respectively
    return (j() < lambda.j() ||
	    (j() == lambda.j() && ((index0().e() < lambda.index0().e() ||
				    (index0().e() == lambda.index0().e() && index1().e() < lambda.index1().e())) ||
				   ((index0().e() == lambda.index0().e() && index1().e() == lambda.index1().e()) &&
				    (index0().k() < lambda.index0().k() ||
				     (index0().k() == lambda.index0().k() && index1().k() < lambda.index1().k()))))));
  }

  template <class BASIS0, class BASIS1>
  TensorProductIndex<BASIS0,BASIS1>&
  TensorProductIndex<BASIS0,BASIS1>::operator ++ ()
  {
    if (index1() == BASIS1::last_index(j(), index1().e())) {
      if (index0() == BASIS0::last_index(j(), index0().e())) {
	if (index1() == BASIS1::last_wavelet(j())) {
 	  if (index0() == BASIS0::last_wavelet(j())) {
 	    index0_ = BASIS0::first_generator(j()+1); // increments j
	    index1_ = BASIS1::first_wavelet(j()); // no generators on higher scales
	  } else {
	    ++index0_;
	    index1_ = BASIS1::first_generator(j());
	  }
	} else {
	  ++index1_;
	  index0_ = BASIS0::first_index(j(), index0().e());
	}
      } else {
	++index0_;
	index1_ = BASIS1::first_index(j(), index1().e());
      }
    } else
      ++index1_;
    
    return *this;
  }

//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   first_generator(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS0,BASIS1>(first_generator(&basis->basis1(), j),
// 					     first_generator(&basis->basis2(), j));
//   }

//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   last_generator(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS0,BASIS1>(last_generator(&basis->basis1(), j),
// 					     last_generator(&basis->basis2(), j));
//   }

//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   first_wavelet(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS0,BASIS1>(basis,
// 					     first_generator(&basis->basis1(), j),
// 					     first_wavelet(&basis->basis2(), j));
//   }

//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   last_wavelet(const int j)
//   {
//     assert(j >= basis->j0());
//     return TensorProductIndex<BASIS0,BASIS1>(last_wavelet(&basis->basis1(), j),
// 					     last_wavelet(&basis->basis2(), j));
//   }
}

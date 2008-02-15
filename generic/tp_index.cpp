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

  template <class BASIS0, class BASIS1>
  const int
  TensorProductIndex<BASIS0,BASIS1>::number() const
  {
    int result = 0;
    
    typedef TensorProductBasis<BASIS0,BASIS1> Basis;

    const int j0 = Basis::j0();

    // determine how many wavelets there are with levels j0<=j'<j
    if (j() > j0)
      result += Basis::Deltasize(j());
   
    // determine how many wavelets there are with level j and type 0<=e'<e
    const int ecode = index0_.e()*2+index1_.e();
    switch(ecode) {
    case 1: // e=(0,1)
      if (j() == j0)
	result +=
	  BASIS0::Deltasize(j())
	  * BASIS1::Deltasize(j());
      break;
    case 2: // e=(1,0)
      result +=
	BASIS0::Deltasize(j())
	* (j() == j0
	   ? BASIS1::Deltasize(j()+1)
	   : BASIS1::Nablasize(j()));
      break;
    case 3: // e=(1,1)
      result += (j() == j0
		 ? BASIS0::Deltasize(j())
		 * BASIS1::Deltasize(j()+1)
		 + BASIS0::Nablasize(j())
		 * BASIS1::Deltasize(j())
		 : BASIS0::Deltasize(j())
		 * BASIS1::Nablasize(j())
		 + BASIS0::Nablasize(j())
		 * BASIS1::Deltasize(j()));
      break;
    case 0: // add nothing
      break;
    }

    // determine how many wavelets there are with level j, type e and k'<k
    switch(ecode) {
    case 0: // e=(0,0)
      result +=
	(index0_.k()-BASIS0::DeltaLmin())
	* BASIS1::Deltasize(j())
	+ index1_.k()-BASIS1::DeltaLmin();
      break;
    case 1: // e=(0,1)
      result +=
	(index0_.k()-BASIS0::DeltaLmin())
	* BASIS1::Nablasize(j())
	+ index1_.k()-BASIS1::Nablamin();
      break;
    case 2:
      result +=
	(index0_.k()-BASIS0::Nablamin())
	* BASIS1::Deltasize(j())
	+ index1_.k()-BASIS1::DeltaLmin();	
      break;
    case 3:
      result +=
	(index0_.k()-BASIS0::Nablamin())
	* BASIS1::Nablasize(j())
	+ index1_.k()-BASIS1::Nablamin();
      break;
    }
    
    return result;
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

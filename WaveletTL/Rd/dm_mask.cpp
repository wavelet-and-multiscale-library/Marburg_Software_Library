// implementation for dm_mask.h

#include <utils/multiindex.h>

namespace WaveletTL
{
  template <class MASK0, class MASK1>
  DMMask1<MASK0, MASK1>::DMMask1()
  {
    MASK0 mask0;
    MASK1 mask1;
    
    for (typename MASK0::const_iterator it0(mask0.begin());
	 it0 != mask0.end(); ++it0)
      for (typename MASK1::const_iterator it1(mask1.begin());
	   it1 != mask1.end(); ++it1)
	{
	  MultivariateLaurentPolynomial<double, 1>::operator []
	    (MultiIndex<int, 1>(it0.index()[0]-it1.index()[0]))
	    += 0.5 * *it0 * *it1;
	}
  }
  
  template <class MASK0, class MASK1, class MASK2>
  DMMask2<MASK0, MASK1, MASK2>::DMMask2()
  {
    MASK0 mask0;
    MASK1 mask1;
    MASK2 mask2;

    for (typename MASK0::const_iterator it0(mask0.begin());
	 it0 != mask0.end(); ++it0)
      for (typename MASK1::const_iterator it1(mask1.begin());
	   it1 != mask1.end(); ++it1)
	for (typename MASK2::const_iterator it2(mask2.begin());
	     it2 != mask2.end(); ++it2)
	  {
	    MultivariateLaurentPolynomial<double, 2>::operator []
	      (MultiIndex<int, 2>(it0.index()[0]-it1.index()[0],
				  it0.index()[0]-it2.index()[0]))
	      += 0.5 * *it0 * *it1 * *it2;
	  }
  }
}

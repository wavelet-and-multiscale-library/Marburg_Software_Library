// implementation for dm_mask.h

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
	  LaurentPolynomial<double>::operator [] (it0.index()-it1.index())
	    += 0.5 * *it0 * *it1;
	}
  }
}

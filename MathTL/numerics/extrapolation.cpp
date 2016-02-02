// implementation for extrapolation.h

#include <cassert>
#include <cmath>

namespace MathTL
{
  template <class ALGORITHM, class RESULT, class SEQUENCE>
  ExtrapolationTable<ALGORITHM,RESULT,SEQUENCE>::ExtrapolationTable
  (const ALGORITHM& T, const unsigned int size, const int p)
    : table_(size, size)
  {
    assert(size >= 2);

    SEQUENCE s;

    // fill first column of the extrapolation table
    for (unsigned int j = 1; j <= size; j++)
      T.approximate(s.n(j), table_(j-1, 0));
    
    // compute right half of the extrapolation table (Aitken-Neville)
    for (unsigned int k = 2; k <= size; k++)
      for (unsigned int j = k; j <= size; j++)
	table_(j-1,k-1) =
	  table_(j-1,k-2)
	  + 1./(pow((double)s.n(j)/(double)s.n(j-k+1),p)-1.) * (table_(j-1,k-2)-table_(j-2,k-2));
  }
}

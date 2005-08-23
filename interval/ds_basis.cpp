// implementation for ds_basis

#include <cassert>

namespace WaveletTL
{
  template <int d, int dT, int ellTl, int ellTr>
  DSBasis<d,dT,ellTl,ellTr>::DSBasis(DSZ bc, const int s, const int sT,
				     DSBiorthogonalizationMethod bio)
  {
    assert(s < d && sT < dT);
    
    switch(bc)
      {
      case empty:
	s0 = s1 = 0;
	sT0 = sT1 = sT;
	break;
      case Zero:
	s0 = s;
	s1 = sT0 = 0;
	sT1 = sT;
	break;
      case One:
	s0 = sT1 = 0;
	s1 = s;
	sT0 = sT;
	break;
      case ZeroOne:
      default:
	s0 = s1 = s;
	sT0 = sT1 = 0;
	break;
      }
  }
}

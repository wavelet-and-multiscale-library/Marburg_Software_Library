// implementation for daubechies_mask.h

#include <cmath>

namespace WaveletTL
{
  template <int N>
  DaubechiesMask<N>::DaubechiesMask()
  {
    const double theta = sqrt(5+2*sqrt(10.));
    switch(N) {
    case 2:
      set_coefficient(MultiIndex<int, 1>(0), (1+sqrt(3.))/4);
      set_coefficient(MultiIndex<int, 1>(1), (3+sqrt(3.))/4);     
      set_coefficient(MultiIndex<int, 1>(2), (3-sqrt(3.))/4);     
      set_coefficient(MultiIndex<int, 1>(3), (1-sqrt(3.))/4);     
      break;
    case 3:
      set_coefficient(MultiIndex<int, 1>(0), (1+sqrt(10.)+theta)/16);
      set_coefficient(MultiIndex<int, 1>(1), (5+sqrt(10.)+3*theta)/16);     
      set_coefficient(MultiIndex<int, 1>(2), (5-sqrt(10.)+theta)/8);     
      set_coefficient(MultiIndex<int, 1>(3), (5-sqrt(10.)-theta)/8);     
      set_coefficient(MultiIndex<int, 1>(4), (5+sqrt(10.)-3*theta)/16);     
      set_coefficient(MultiIndex<int, 1>(5), (1+sqrt(10.)-theta)/16);     
      break;
    case 1:
    default: // Haar mask
      set_coefficient(MultiIndex<int, 1>(0), 1);
      set_coefficient(MultiIndex<int, 1>(1), 1);
    }
  }

}

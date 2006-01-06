// implementation for p_basis.h

#include <cassert>

namespace WaveletTL
{
  template <int d, int dT>
  PBasis<d,dT>::PBasis(const int s0, const int s1, const int sT0, const int sT1) {
    assert(std::max(s0,s1) < d && std::max(sT0,sT1) < dT);
        
    this->s0 = s0;
    this->s1 = s1;
    this->sT0 = sT0;
    this->sT1 = sT1;

    setup();
  }


  template <int d, int dT>
  void
  PBasis<d,dT>::setup() {
    j0_ = (int) ceil(log((double)d)/M_LN2); // for the moment
  }

}

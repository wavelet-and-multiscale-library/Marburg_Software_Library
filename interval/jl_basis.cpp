// implementation for jl_basis.h

namespace WaveletTL
{
  JLBasis::JLBasis(const int s0, const int s1)
  {
    this->s0 = s0;
    this->s1 = s1;
    
    setup();
  }
  
  void
  JLBasis::setup() {
    j0_ = 1;
  }
  
}

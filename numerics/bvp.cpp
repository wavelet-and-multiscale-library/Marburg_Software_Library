// implementation for bvp.h

namespace MathTL
{
  template <unsigned int DIM>
  inline
  TwoPointBVP<DIM>::~TwoPointBVP()
  {
  }

  template <class ATLAS>
  EllipticBVP<ATLAS>::EllipticBVP(const ATLAS& atlas)
    : atlas_(atlas)
  {
  }
}

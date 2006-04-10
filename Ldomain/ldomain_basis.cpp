// implementation for ldomain_basis.h

namespace WaveletTL
{
  template <class IBASIS>
  LDomainBasis<IBASIS>::LDomainBasis()
    : basis01(true, false),
      basis10(false, true)
  {
  }

}

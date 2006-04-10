// implementation for ldomain_basis.h

namespace WaveletTL
{
  template <class IBASIS>
  LDomainBasis<IBASIS>::LDomainBasis()
    : basis00_(true, true),
      basis01_(true, false),
      basis10_(false, true)
  {
  }

}

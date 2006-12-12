// implementation for ldomain_spline_basis.h

namespace WaveletTL
{
  template <int d, int dT>
  LDomainSplineBasis<d,dT>::LDomainSplineBasis()
    : basis1D_01("P","",0,1,0,0),
      basis1D_10("P","",1,0,0,0)
  {
  }
}

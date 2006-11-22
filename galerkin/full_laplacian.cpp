// implementation for full_laplacian.h

namespace WaveletTL
{
  template <int d, int dT>
  FullLaplacian<d,dT>::FullLaplacian(const SplineBasisData<d,dT>& sd)
    : sd_(sd)
  {
  }

}

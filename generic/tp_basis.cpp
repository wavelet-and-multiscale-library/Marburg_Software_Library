// implementation for tp_basis.h

namespace WaveletTL
{
  template <class BASIS1, class BASIS2>
  TensorProductBasis<BASIS1,BASIS2>::TensorProductBasis()
  {
    j0_ = std::max(basis1_.j0(), basis2_.j0());
  }
}

// implementation for tp_index.h

namespace WaveletTL
{
  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>::TensorProductIndex(const BASIS1* basis1,
							const BASIS2* basis2)
    : basis1_(basis1), basis2_(basis2), index1_(basis1), index2_(basis2)
  {
  }

}

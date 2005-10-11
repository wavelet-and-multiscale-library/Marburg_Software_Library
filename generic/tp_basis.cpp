// implementation for tp_basis.h

namespace WaveletTL
{
  template <class BASIS1, class BASIS2>
  TensorProductBasis<BASIS1,BASIS2>::TensorProductBasis()
  {
    j0_ = std::max(basis1_.j0(), basis2_.j0());
  }

  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::decompose(const InfiniteVector<double, Index>& c,
					       const int j0,
					       InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_1(it.index(), j0, help);
      d += *it * help;
    }
  }
  
  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::decompose_t(const InfiniteVector<double, Index>& c,
						 const int j0,
						 InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      decompose_t_1(it.index(), j0, help);
      d += *it * help;
    }
  }

  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::reconstruct(const InfiniteVector<double, Index>& c,
						 const int j,
						 InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d += *it * help;
    }
  }
  
  template <class BASIS1, class BASIS2>
  void
  TensorProductBasis<BASIS1,BASIS2>::reconstruct_t(const InfiniteVector<double, Index>& c,
						   const int j,
						   InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      d += *it * help;
    }
  }

}

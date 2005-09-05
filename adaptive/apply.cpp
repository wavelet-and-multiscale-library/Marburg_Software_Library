// implementation for APPLY

namespace WaveletTL
{
  template <class PROBLEM>
  void APPLY(const PROBLEM& P,
	     const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& v,
	     const double eta,
	     InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& w)
  {
    w.clear();

    typedef InfiniteVector<double, typename PROBLEM::WaveletBasis::Index> VectorType;
    for (typename VectorType::const_iterator it(v.begin()); it != v.end(); ++it) {
      // quick hack, only to make add_column compile
      P.add_column(*it, it.index(), 2, w);
    }
  }  
}

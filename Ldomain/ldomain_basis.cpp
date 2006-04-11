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

  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct(const InfiniteVector<double, Index>& c,
				    const int j,
				    InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_t(const InfiniteVector<double, Index>& c,
				      const int j,
				      InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      d.add(*it, help);
    }
  }

  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_1(const Index& lambda,
				      const int j,
				      InfiniteVector<double, Index>& c) const {
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_t_1(const Index& lambda,
					const int j,
					InfiniteVector<double, Index>& c) const {
  }

}

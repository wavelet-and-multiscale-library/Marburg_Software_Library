// implementation for ldomain_gramian.h

namespace WaveletTL
{
  template <int d, int dT>
  LDomainGramian<SplineBasis<d,dT,DS_construction> >::LDomainGramian
  (const WaveletBasis& basis,
   const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y), normA(0.0), normAinv(0.0)
  {
  }

}

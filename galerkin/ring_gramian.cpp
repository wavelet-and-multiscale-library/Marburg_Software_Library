// implementation for ring_gramian.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingGramian<d,dt,s0,s1>::RingGramian(const WaveletBasis& basis,
				       const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y)
  {
  }


}

// implementation for ring_gramian.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingGramian<d,dt,s0,s1>::RingGramian(const RingBasis<d,dt,s0,s1>& basis,
				       const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y)
  {
  }
  
  template <int d, int dt, int s0, int s1>
  inline
  double
  RingGramian<d,dt,s0,s1>::a(const Index& lambda,
			     const Index& nu) const
  {
    return basis_.integrate(lambda, nu);
  }

}

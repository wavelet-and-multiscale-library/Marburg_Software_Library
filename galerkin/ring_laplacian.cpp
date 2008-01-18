// implementation for ring_laplacian.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingLaplacian<d,dt,s0,s1>::RingLaplacian(const RingBasis<d,dt,s0,s1>& basis,
					   const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y)
  {
  }
  
  template <int d, int dt, int s0, int s1>
  inline
  double
  RingLaplacian<d,dt,s0,s1>::a(const Index& lambda,
			       const Index& nu) const
  {
    // just to let the test program compile
    return basis_.integrate(lambda, nu);
  }
  
}

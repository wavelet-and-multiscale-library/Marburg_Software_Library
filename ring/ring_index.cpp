// implementation for ring.h

#include <interval/periodic.h>
#include <interval/spline_basis.h>

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingIndex<d,dt,s0,s1>::RingIndex()
    : j_(RingBasis<d,dt,s0,s1>::j0()),
      k_(PeriodicBasis<CDFBasis<d,dt> >::DeltaLmin(),
	 SplineBasis<d,dt,P_construction,s0,s1,0,0>::DeltaLmin())
  {
  }
}

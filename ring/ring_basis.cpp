// implementation for ring_basis.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::first_generator(const int j)
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::DeltaLmin(),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::DeltaLmin());
    
    return Index(j, e, k);
  }
  
  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::last_generator(const int j)
  {
    assert(j >= j0());

    typename Index::type_type e;
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::DeltaRmax(j),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::DeltaRmax(j));
    
    return Index(j, e, k);
  }

  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::first_wavelet(const int j)
  {
    assert(j >= j0());
    
    typename Index::type_type e(0, 1);
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::DeltaLmin(),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::Nablamin());

    return Index(j, e, k);
  }

  template <int d, int dt, int s0, int s1>
  typename RingBasis<d,dt,s0,s1>::Index
  RingBasis<d,dt,s0,s1>::last_wavelet(const int j)
  {
    assert(j >= j0());
    
    typename Index::type_type e(1, 1);
    typename Index::translation_type k
      (PeriodicBasis<CDFBasis<d,dt> >::Nablamax(j),
       SplineBasis<d,dt,P_construction,s0,s1,0,0>::Nablamax(j));
    
    return Index(j, e, k);
  }
  
  
}

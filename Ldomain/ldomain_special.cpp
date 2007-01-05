// implementation for ldomain_basis.h, template specialization to SplineBasis<d,dT,DS_construction>

#include <cmath>
#include <time.h>
#include <iostream>

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <int d, int dT>
  LDomainBasis<SplineBasis<d,dT,DS_construction> >
  ::LDomainBasis(const IntervalBasis& basis1d)
    : basis1d_(basis1d)
  {
// #if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
//     supp_hits = 0;
//     supp_misses = 0;
//     Mj1_hits = 0;
//     Mj1_misses = 0;
// #endif
  }

  template <int d, int dT>
  const int
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::Deltasize(const int j) const {
    const unsigned int Deltaj = basis1d().Deltasize(j);
    return 3*(Deltaj-2)*(Deltaj-2)+2*(Deltaj-2);
  }
  
  template <int d, int dT>
  const int
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::Nabla01size(const int j) const {
    return (3*basis1d().Deltasize(j)-5)*(1<<j);
  }
  
  template <int d, int dT>
  const int
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::Nabla10size(const int j) const {
    return (3*basis1d().Deltasize(j)-5)*(1<<j);
  }
  
  template <int d, int dT>
  const int
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::Nabla11size(const int j) const {
    return 3*(1<<(2*j));
  }
  
  template <int d, int dT>
  typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::first_generator(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e;

    // setup lowest translation index for e=(0,0), p=0
    typename Index::translation_type k(basis1d().DeltaLmin()+1, basis1d().DeltaLmin()+1);
    
    return Index(j, e, 0, k, this);
  }

  template <int d, int dT>
  typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::last_generator(const int j) const
  {
    assert(j >= j0());
    
    typename Index::type_type e;
    
    // setup highest translation index for e=(0,0), p=4
    typename Index::translation_type k(0, basis1d().DeltaRmax(j)-1);
    
    return Index(j, e, 4, k, this);
  }

  template <int d, int dT>
  typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::first_wavelet(const int j) const
  {
    assert(j >= j0());

    typename Index::type_type e(0, 1);

    // setup lowest translation index for e=(0,1), p=0
    typename Index::translation_type k(basis1d().DeltaLmin()+1, basis1d().Nablamin());
    
    return Index(j, e, 0, k, this);
  }

  template <int d, int dT>
  typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index
  LDomainBasis<SplineBasis<d,dT,DS_construction> >
  ::first_wavelet(const int j, const typename Index::type_type& ewish) const
  {
    assert(j >= j0());
    
    typename Index::type_type e(ewish);
    
    // setup lowest translation index appropriately
    typename Index::translation_type k;
    const int ecode(e[0]+2*e[1]);
    if (ecode == 0) {
      // e = (0,0)
      k[0] = k[1] = basis1d().DeltaLmin()+1;
    } else {
      if (ecode == 1) {
	// e = (1,0)
	k[0] = basis1d().Nablamin();
	k[1] = basis1d().DeltaLmin()+1;
      } else {
	if (ecode == 2) {
	  // e = (0,1)
	  k[0] = basis1d().DeltaLmin()+1;
	  k[1] = basis1d().Nablamin();
	} else {
	  // e = (1,1)
	  k[0] = k[1] = basis1d().Nablamin();
	}
      }
    }
    
    return Index(j, e, 0, k, this);
  }

  template <int d, int dT>
  typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index  
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::last_wavelet(const int j) const
  {
    assert(j >= j0());
    
    typename Index::type_type e(1, 1);
    
    // setup highest translation index for e=(1,1), p=2
    typename Index::translation_type k(basis1d().Nablamax(j), basis1d().Nablamax(j));
    
    return Index(j, e, 2, k, this);
  }


}

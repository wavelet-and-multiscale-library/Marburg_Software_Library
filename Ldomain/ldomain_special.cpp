// implementation for ldomain_basis.h, template specialization to SplineBasis<d,dT,DS_construction>

#include <cmath>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <utils/map_tools.h>

using std::cout;
using std::endl;

namespace WaveletTL
{
  template <int d, int dT>
  LDomainBasis<SplineBasis<d,dT,DS_construction> >
  ::LDomainBasis(const IntervalBasis& basis1d)
    : basis1d_(basis1d)
  {
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
    supp_hits = 0;
    supp_misses = 0;
    rec1_hits = 0;
    rec1_misses = 0;
#endif
  }

  template <int d, int dT>
  const int
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::Deltasize(const int j) const {
    const unsigned int Deltaj = basis1d().Deltasize(j);
    return (3*(Deltaj-2)+2)*(Deltaj-2);
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

  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::support
  (const Index& lambda, Support& supp) const
  {
    // check whether the supp(psi_lambda) already exists in the cache
    typename SupportCache::iterator supp_lb(supp_cache.lower_bound(lambda));
    typename SupportCache::iterator supp_it(supp_lb);
    if (supp_lb == supp_cache.end() ||
	supp_cache.key_comp()(lambda, supp_lb->first))
      {
	// compute supp(psi_lambda) and insert it into the cache
	typedef typename SupportCache::value_type value_type;

	const int ecode = lambda.e()[0]+2*lambda.e()[1];
	const int lambdaj = lambda.j();
	
	if (ecode == 0) {
	  // psi_lambda is a generator. Here we know by construction of the
	  // composite basis that per patch, psi_lambda looks like a single
	  // tensor product of 1D generators (possibly weighted by a factor).
	  
	  supp.j = lambdaj;
	  
	  switch (lambda.p()) {
	  case 0:
	    // psi_lambda completely lives on patch 0
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[0],
							    &basis1d()),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[1],
							    &basis1d()),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    supp.xmin[1] = supp.xmin[2] = -1;
	    
	    break;
	  case 1:
	    // psi_lambda completely lives on patch 1
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[0],
							    &basis1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[1],
							    &basis1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[0] = supp.xmin[2] = -1;
	    
	    break;
	  case 2:
	    // psi_lambda completely lives on patch 2
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj, 
							    0,
							    lambda.k()[0],
							    &basis1d()),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0, 
							    lambda.k()[1],
							    &basis1d()),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = supp.xmin[1] = -1;
	    
	    break;
	  case 3:
	    // psi_lambda lives on patches 0 and 1
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[0],
							    &basis1d()),
			      supp.xmin[0],
			      supp.xmax[0]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    basis1d().DeltaLmin(),
							    &basis1d()),
			      supp.ymin[0],
			      supp.ymax[0]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[0],
							    &basis1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    basis1d().DeltaRmax(lambdaj),
							    &basis1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    supp.xmin[2] = -1;
	    
	    break;
	  case 4:
	    // psi_lambda lives on patches 1 and 2
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    basis1d().DeltaRmax(lambdaj),
							    &basis1d()),
			      supp.xmin[1],
			      supp.xmax[1]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[1],
							    &basis1d()),
			      supp.ymin[1],
			      supp.ymax[1]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    basis1d().DeltaLmin(),
							    &basis1d()),
			      supp.xmin[2],
			      supp.xmax[2]);
	    
	    basis1d().support(typename IntervalBasis::Index(lambdaj,
							    0,
							    lambda.k()[1],
							    &basis1d()),
			      supp.ymin[2],
			      supp.ymax[2]);
	    
	    supp.xmin[0] = -1;
	    
	    break;
	  }
	} else {
	  // wavelet
	  
	  supp.j = lambdaj+1;
	  
	  // compute the expansion coefficients of psi_lambda w.r.t. the
	  // generators of the next higher scale, then aggregating all the supports
	  // (of course, this is a brute force solution...)
	  
 	  InfiniteVector<double, Index> gcoeffs; // dummy object to let the following line compile
	  typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()), gcoeffs_end(gcoeffs.begin());
	  reconstruct_1(lambda, lambdaj+1, it, gcoeffs_end);
	  
	  Support tempsupp;
	  
	  // initialize the support with an "empty" set
	  for (int p = 0; p <= 2; p++) {
	    supp.xmin[p] = -1;
	  }
	  
	  for (; it != gcoeffs_end; ++it)
	    {
	      // compute supp(psi_mu)
	      support(it.index(), tempsupp);

	      // unite the two support sets
	      for (int p = 0; p <= 2; p++) {
		if (tempsupp.xmin[p] != -1) {
		  // a nontrivial new support share, we have to do something
		  if (supp.xmin[p] == -1) {
		    // previous support estimate was "empty", we have to insert a nontrivial new one
		    supp.xmin[p] = tempsupp.xmin[p];
		    supp.xmax[p] = tempsupp.xmax[p];
		    supp.ymin[p] = tempsupp.ymin[p];
		    supp.ymax[p] = tempsupp.ymax[p];
		  } else {
		    // previous support estimate was nontrivial, we have to compute a new one
		    supp.xmin[p] = std::min(tempsupp.xmin[p], supp.xmin[p]);
		    supp.xmax[p] = std::max(tempsupp.xmax[p], supp.xmax[p]);
		    supp.ymin[p] = std::min(tempsupp.ymin[p], supp.ymin[p]);
		    supp.ymax[p] = std::max(tempsupp.ymax[p], supp.ymax[p]);
		  }
		}
	      }     
	    }
	}
	
	supp_it = supp_cache.insert(supp_lb, value_type(lambda, supp));

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	supp_misses++;
	if ((supp_hits+supp_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis support cache (hits/misses/total/hit ratio): ("
		 << supp_hits << "/"
		 << supp_misses << "/"
		 << supp_hits+supp_misses << "/"
		 << (double) supp_hits/(supp_hits+supp_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif
      }
    else
      {
	// cache hit, copy the precomputed support
	const Support& suppcache = supp_it->second;
  	supp.j = suppcache.j;
  	for (unsigned int i = 0; i < 3; i++) {
  	  supp.xmin[i] = suppcache.xmin[i];
  	  supp.xmax[i] = suppcache.xmax[i];
  	  supp.ymin[i] = suppcache.ymin[i];
  	  supp.ymax[i] = suppcache.ymax[i];
  	}
	
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	supp_hits++;
	if ((supp_hits+supp_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis support cache (hits/misses/total/hit ratio): ("
		 << supp_hits << "/"
		 << supp_misses << "/"
		 << supp_hits+supp_misses << "/"
		 << (double) supp_hits/(supp_hits+supp_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif
      }  
  }
  
  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::reconstruct_1
  (const Index& lambda,
   const int j,
   InfiniteVector<double, Index>& c) const
  {
    assert(j == lambda.j()+1); // we only support a level difference of one
    
    typedef std::map<size_type,double> V;

    // check whether the result already exists in the cache
    typename Reconstruct1Cache::iterator rec1_lb(rec1_cache.lower_bound(lambda));
    typename Reconstruct1Cache::iterator rec1_it(rec1_lb);
    if (rec1_lb == rec1_cache.end() ||
	rec1_cache.key_comp()(lambda, rec1_lb->first))
      {
     	// compute supp(psi_lambda) and insert it into the cache
	typedef typename Reconstruct1Cache::value_type value_type;

	const unsigned int Deltaj   = basis1d().Deltasize(j);
	const MultiIndex<int,2> zero(0,0); // e=(0,0)
	
	const int ecode(lambda.e()[0]+2*lambda.e()[1]);
	switch(ecode) {
	case 0:
	  // generator
	  basis1d().Mj0_.set_level(lambda.j());
	  switch(lambda.p()) {
	  case 0: {
	    // psi_lambda decomposes into generators on patch 0 alone
	    // apply kron(M#,M#), cf. KroneckerMatrix::apply()
	    // 1. first decide which block number and subindex corresponds to the given generator
	    const size_type block_nr    = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    
	    // 2. get column of second factor M#
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    // 3. get column of first factor M#
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    // 4. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					0,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+1+it2->first),
					this),
				  it2->second * it3->second);	
	  }
	    break;
	  case 1: {
	    // psi_lambda decomposes into generators on patch 1 alone
	    // apply kron(M#,M#)
	    const size_type block_nr    = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					1,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+1+it2->first),
					this),
				  it2->second * it3->second);	
	  }
	    break;
	  case 2: {
	    // psi_lambda decomposes into generators on patch 2 alone
	    // apply kron(M#,M#)
	    const size_type block_nr    = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					2,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+1+it2->first),
					this),
				  it2->second * it3->second);	
	  }
	    break;
	  case 3: {
	    // psi_lambda decomposes into generators on patches 0,1 and 3
	    
	    // contribution from patch 0
	    //
	    // apply kron(M#,ML):
	    // 1. get column of second factor ML
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(0, 1.0));
	    basis1d().Mj0_.apply(z1, z2); // we have to neglect the first entry of z2 later
	    
	    // 2. get column of first factor M#
	    z1.clear();
	    const size_type block_nr = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    // 3. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(++(z2.begin())); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					0,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+it2->first), // note the "-1"
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 1
	    //
	    // apply kron(M#,MR):
	    // 1. get column of second factor MR
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
	    z2.clear();
	    basis1d().Mj0_.apply(z1, z2); // we have to neglect the last entry of z2 later
	    
	    // 2. get column of first factor M#: this has already been done above
	    
	    // 3. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::reverse_iterator it2(++(z2.rbegin())); it2 != z2.rend(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					1,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+it2->first), // note the "-1"
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 3
	    //
	    // apply kron(M#,Mtopleft):
	    // 1. get column of first factor M#: this has already been done above
	    
	    // 2. get top left entry of Mj0
	    const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	    
	    // 3. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      c.add_coefficient(Index(lambda.j()+1,
				      zero,
				      3,
				      MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							0),
				      this),
				it3->second * Mtopleft);
	  }
	    break;
	  case 4: {
	    // psi_lambda decomposes into generators on patches 1,2 and 4
	    
	    // contribution from patch 1
	    //
	    // apply kron(MR,M#):
	    // 1. get column of second factor M#
	    V z1, z2, z3;
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    // 2. get column of first factor MR
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
	    basis1d().Mj0_.apply(z1, z3); // we have to neglect the last entry of z3 later
	    
	    // 3. combine results
	    for (V::reverse_iterator it3(++(z3.rbegin())); it3 != z3.rend(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					1,
					MultiIndex<int,2>(basis1d().DeltaLmin()+it3->first, // note the "-1"
							  basis1d().DeltaLmin()+1+it2->first), 
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 2
	    //
	    // apply kron(ML,M#):
	    // 1. get column of second factor M#: this has already been done above
	    
	    // 2. get column of first factor ML
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(0, 1.0));
	    z3.clear();
	    basis1d().Mj0_.apply(z1, z3); // we have to neglect the first entry of z3 later
	    
	    // 3. combine results
	    for (V::const_iterator it3(++(z3.begin())); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					2,
					MultiIndex<int,2>(basis1d().DeltaLmin()+it3->first, // note the "-1"
							  basis1d().DeltaLmin()+1+it2->first), 
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 4
	    //
	    // apply kron(Mtopleft,M#):
	    // 1. get column of second factor M#: this has already been done above
	    
	    // 2. get top left entry of Mj0
	    const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	    
	    // 3. combine results
	    for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	      c.add_coefficient(Index(lambda.j()+1,
				      zero,
				      4,
				      MultiIndex<int,2>(0,
							basis1d().DeltaLmin()+1+it2->first),
				      this),
				it2->second * Mtopleft);
	  }
	    break;
	  default:
	    break;
	  }
	  break;
	case 1:
	  // (1,0)-wavelet
	  {
	    std::map<size_type,double> x,y,z,z1,z2;
	    x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j())-Nabla01size(lambda.j()), 1.0));
	    apply_Mj1c_10(lambda.j(), x, y); // generator coeffs of initial stable completion
	    apply_Mj0T_transposed(lambda.j(),  y, z1);
	    apply_Mj0            (lambda.j(), z1, z2);
	    add_maps(y, z2, z, 1.0, -1.0);
	    map_to_vector(lambda.j()+1, z, c);
	  }
	  break;
	case 2:
	  // (0,1)-wavelet
	  {
	    std::map<size_type,double> x,y,z,z1,z2;
	    x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j()), 1.0));
	    apply_Mj1c_01(lambda.j(), x, y); // generator coeffs of initial stable completion
	    apply_Mj0T_transposed(lambda.j(),  y, z1);
	    apply_Mj0            (lambda.j(), z1, z2);
	    add_maps(y, z2, z, 1.0, -1.0);
	    map_to_vector(lambda.j()+1, z, c);
	  }
	  break;
	case 3:
	  // (1,1)-wavelet
	  {
	    std::map<size_type,double> x,y,z,z1,z2;
	    x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j())-Nabla01size(lambda.j())-Nabla10size(lambda.j()), 1.0));
	    apply_Mj1c_11(lambda.j(), x, y); // generator coeffs of initial stable completion
	    apply_Mj0T_transposed(lambda.j(),  y, z1);
	    apply_Mj0            (lambda.j(), z1, z2);
	    add_maps(y, z2, z, 1.0, -1.0);
	    map_to_vector(lambda.j()+1, z, c);
	  }
	  break;
	default:
	  break;
	}

	rec1_it = rec1_cache.insert(rec1_lb, value_type(lambda, c));
	
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	rec1_misses++;
	if ((rec1_hits+rec1_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis reconstruct_1 cache (hits/misses/total/hit ratio): ("
		 << rec1_hits << "/"
		 << rec1_misses << "/"
		 << rec1_hits+rec1_misses << "/"
		 << (double) rec1_hits/(rec1_hits+rec1_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif

      }
    else
      {
	// cache hit, copy the precomputed coefficients
	c = rec1_it->second;

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	rec1_hits++;
	if ((rec1_hits+rec1_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis reconstruct_1 cache (hits/misses/total/hit ratio): ("
		 << rec1_hits << "/"
		 << rec1_misses << "/"
		 << rec1_hits+rec1_misses << "/"
		 << (double) rec1_hits/(rec1_hits+rec1_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif
      }
  }

  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::reconstruct_1
  (const Index& lambda, const int j,
   typename InfiniteVector<double,Index>::const_iterator& it_begin,
   typename InfiniteVector<double,Index>::const_iterator& it_end) const
  {
    assert(j == lambda.j()+1); // we only support a level difference of one
    
    typedef std::map<size_type,double> V;
    
    // check whether the result already exists in the cache
    typename Reconstruct1Cache::iterator rec1_lb(rec1_cache.lower_bound(lambda));
    typename Reconstruct1Cache::iterator rec1_it(rec1_lb);
    if (rec1_lb == rec1_cache.end() ||
	rec1_cache.key_comp()(lambda, rec1_lb->first))
      {
     	// compute supp(psi_lambda) and insert it into the cache
	typedef typename Reconstruct1Cache::value_type value_type;

	InfiniteVector<double,Index> c;

	const unsigned int Deltaj   = basis1d().Deltasize(j);
	const MultiIndex<int,2> zero(0,0); // e=(0,0)
	
	const int ecode(lambda.e()[0]+2*lambda.e()[1]);
	switch(ecode) {
	case 0:
	  // generator
	  basis1d().Mj0_.set_level(lambda.j());
	  switch(lambda.p()) {
	  case 0: {
	    // psi_lambda decomposes into generators on patch 0 alone
	    // apply kron(M#,M#), cf. KroneckerMatrix::apply()
	    // 1. first decide which block number and subindex corresponds to the given generator
	    const size_type block_nr    = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    
	    // 2. get column of second factor M#
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    // 3. get column of first factor M#
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    // 4. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					0,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+1+it2->first),
					this),
				  it2->second * it3->second);	
	  }
	    break;
	  case 1: {
	    // psi_lambda decomposes into generators on patch 1 alone
	    // apply kron(M#,M#)
	    const size_type block_nr    = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					1,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+1+it2->first),
					this),
				  it2->second * it3->second);	
	  }
	    break;
	  case 2: {
	    // psi_lambda decomposes into generators on patch 2 alone
	    // apply kron(M#,M#)
	    const size_type block_nr    = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					2,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+1+it2->first),
					this),
				  it2->second * it3->second);	
	  }
	    break;
	  case 3: {
	    // psi_lambda decomposes into generators on patches 0,1 and 3
	    
	    // contribution from patch 0
	    //
	    // apply kron(M#,ML):
	    // 1. get column of second factor ML
	    V z1, z2, z3;
	    z1.insert(std::pair<size_type,double>(0, 1.0));
	    basis1d().Mj0_.apply(z1, z2); // we have to neglect the first entry of z2 later
	    
	    // 2. get column of first factor M#
	    z1.clear();
	    const size_type block_nr = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	    z1.insert(std::pair<size_type,double>(block_nr, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z3);
	    
	    // 3. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(++(z2.begin())); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					0,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+it2->first), // note the "-1"
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 1
	    //
	    // apply kron(M#,MR):
	    // 1. get column of second factor MR
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
	    z2.clear();
	    basis1d().Mj0_.apply(z1, z2); // we have to neglect the last entry of z2 later
	    
	    // 2. get column of first factor M#: this has already been done above
	    
	    // 3. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      for (V::reverse_iterator it2(++(z2.rbegin())); it2 != z2.rend(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					1,
					MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							  basis1d().DeltaLmin()+it2->first), // note the "-1"
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 3
	    //
	    // apply kron(M#,Mtopleft):
	    // 1. get column of first factor M#: this has already been done above
	    
	    // 2. get top left entry of Mj0
	    const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	    
	    // 3. combine results
	    for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	      c.add_coefficient(Index(lambda.j()+1,
				      zero,
				      3,
				      MultiIndex<int,2>(basis1d().DeltaLmin()+1+it3->first,
							0),
				      this),
				it3->second * Mtopleft);
	  }
	    break;
	  case 4: {
	    // psi_lambda decomposes into generators on patches 1,2 and 4
	    
	    // contribution from patch 1
	    //
	    // apply kron(MR,M#):
	    // 1. get column of second factor M#
	    V z1, z2, z3;
	    const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
	    z1.insert(std::pair<size_type,double>(block_index, 1.0));
	    basis1d().Mj0_.apply_central_block(z1, z2);
	    
	    // 2. get column of first factor MR
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(Deltaj-1, 1.0));
	    basis1d().Mj0_.apply(z1, z3); // we have to neglect the last entry of z3 later
	    
	    // 3. combine results
	    for (V::reverse_iterator it3(++(z3.rbegin())); it3 != z3.rend(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					1,
					MultiIndex<int,2>(basis1d().DeltaLmin()+it3->first, // note the "-1"
							  basis1d().DeltaLmin()+1+it2->first), 
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 2
	    //
	    // apply kron(ML,M#):
	    // 1. get column of second factor M#: this has already been done above
	    
	    // 2. get column of first factor ML
	    z1.clear();
	    z1.insert(std::pair<size_type,double>(0, 1.0));
	    z3.clear();
	    basis1d().Mj0_.apply(z1, z3); // we have to neglect the first entry of z3 later
	    
	    // 3. combine results
	    for (V::const_iterator it3(++(z3.begin())); it3 != z3.end(); ++it3)
	      for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
		c.add_coefficient(Index(lambda.j()+1,
					zero,
					2,
					MultiIndex<int,2>(basis1d().DeltaLmin()+it3->first, // note the "-1"
							  basis1d().DeltaLmin()+1+it2->first), 
					this),
				  it2->second * it3->second * M_SQRT1_2);	
	    
	    // contribution from patch 4
	    //
	    // apply kron(Mtopleft,M#):
	    // 1. get column of second factor M#: this has already been done above
	    
	    // 2. get top left entry of Mj0
	    const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	    
	    // 3. combine results
	    for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	      c.add_coefficient(Index(lambda.j()+1,
				      zero,
				      4,
				      MultiIndex<int,2>(0,
							basis1d().DeltaLmin()+1+it2->first),
				      this),
				it2->second * Mtopleft);
	  }
	    break;
	  default:
	    break;
	  }
	  break;
	case 1:
	  // (1,0)-wavelet
	  {
	    std::map<size_type,double> x,y,z,z1,z2;
	    x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j())-Nabla01size(lambda.j()), 1.0));
	    apply_Mj1c_10(lambda.j(), x, y); // generator coeffs of initial stable completion
	    apply_Mj0T_transposed(lambda.j(),  y, z1);
	    apply_Mj0            (lambda.j(), z1, z2);
	    add_maps(y, z2, z, 1.0, -1.0);
	    map_to_vector(lambda.j()+1, z, c);
	  }
	  break;
	case 2:
	  // (0,1)-wavelet
	  {
	    std::map<size_type,double> x,y,z,z1,z2;
	    x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j()), 1.0));
	    apply_Mj1c_01(lambda.j(), x, y); // generator coeffs of initial stable completion
	    apply_Mj0T_transposed(lambda.j(),  y, z1);
	    apply_Mj0            (lambda.j(), z1, z2);
	    add_maps(y, z2, z, 1.0, -1.0);
	    map_to_vector(lambda.j()+1, z, c);
	  }
	  break;
	case 3:
	  // (1,1)-wavelet
	  {
	    std::map<size_type,double> x,y,z,z1,z2;
	    x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j())-Nabla01size(lambda.j())-Nabla10size(lambda.j()), 1.0));
	    apply_Mj1c_11(lambda.j(), x, y); // generator coeffs of initial stable completion
	    apply_Mj0T_transposed(lambda.j(),  y, z1);
	    apply_Mj0            (lambda.j(), z1, z2);
	    add_maps(y, z2, z, 1.0, -1.0);
	    map_to_vector(lambda.j()+1, z, c);
	  }
	  break;
	default:
	  break;
	}

	rec1_it = rec1_cache.insert(rec1_lb, value_type(lambda, c));
	
#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	rec1_misses++;
	if ((rec1_hits+rec1_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis reconstruct_1 cache (hits/misses/total/hit ratio): ("
		 << rec1_hits << "/"
		 << rec1_misses << "/"
		 << rec1_hits+rec1_misses << "/"
		 << (double) rec1_hits/(rec1_hits+rec1_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif

      }
    else
      {
	// cache hit

#if _WAVELETTL_LDOMAINBASIS_VERBOSITY >= 1
	rec1_hits++;
	if ((rec1_hits+rec1_misses)%1000000 == 0)
	  {
	    cout << "[LdomainBasis reconstruct_1 cache (hits/misses/total/hit ratio): ("
		 << rec1_hits << "/"
		 << rec1_misses << "/"
		 << rec1_hits+rec1_misses << "/"
		 << (double) rec1_hits/(rec1_hits+rec1_misses)
		 << ")"
		 << "]"
		 << endl;
	  }
#endif
      }

    it_begin = rec1_it->second.begin();
    it_end = rec1_it->second.end();
  }
  
  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::reconstruct_1
  (const Index& lambda, const int j,
   std::map<size_type,double>& z) const
  {
    assert(j == lambda.j()+1); // we only support a level difference of one
    
    typedef std::map<size_type,double> V;

    const unsigned int Deltaj   = basis1d().Deltasize(j);
    const MultiIndex<int,2> zero(0,0); // e=(0,0)
    
    const int ecode(lambda.e()[0]+2*lambda.e()[1]);
    switch(ecode) {
    case 1:
      // (1,0)-wavelet
      {
	std::map<size_type,double> x,y,z1,z2;
	x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j())-Nabla01size(lambda.j()), 1.0));
	apply_Mj1c_10(lambda.j(), x, y); // generator coeffs of initial stable completion
	apply_Mj0T_transposed(lambda.j(),  y, z1);
	apply_Mj0            (lambda.j(), z1, z2);
 	add_maps(y, z2, z, 1.0, -1.0);
      }
      break;
    case 2:
      // (0,1)-wavelet
      {
	std::map<size_type,double> x,y,z1,z2;
	x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j()), 1.0));
	apply_Mj1c_01(lambda.j(), x, y); // generator coeffs of initial stable completion
	apply_Mj0T_transposed(lambda.j(),  y, z1);
	apply_Mj0            (lambda.j(), z1, z2);
 	add_maps(y, z2, z, 1.0, -1.0);
      }
      break;
    case 3:
      // (1,1)-wavelet
      {
	std::map<size_type,double> x,y,z1,z2;
	x.insert(std::pair<size_type,double>(lambda.number()-Deltasize(lambda.j())-Nabla01size(lambda.j())-Nabla10size(lambda.j()), 1.0));
	apply_Mj1c_11(lambda.j(), x, y); // generator coeffs of initial stable completion
	apply_Mj0T_transposed(lambda.j(),  y, z1);
	apply_Mj0            (lambda.j(), z1, z2);
 	add_maps(y, z2, z, 1.0, -1.0);
      }
      break;
    default: // also for case 0, which we assume not to have to treat here
      break;
    }
  }

  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::map_to_vector
  (const int j,
   const std::map<size_type,double>& x,
   InfiniteVector<double, Index>& c) const
  {
    typedef std::map<size_type,double> V;

    const unsigned int Deltaj   = basis1d().Deltasize(j);
    
    c.clear();
    const MultiIndex<int,2> zero(0,0); // e=(0,0)

    for (V::const_iterator it(x.begin()); it != x.end(); ++it) {
      if (fabs(it->second) > 1e-12) {
	size_type help(it->first);
	unsigned int patch = 0;
	if (help < 3*(Deltaj-2)*(Deltaj-2)) {
	  patch = help / ((Deltaj-2)*(Deltaj-2));
	  help -= patch * (Deltaj-2)*(Deltaj-2);

	  c.add_coefficient(Index(j,
				  zero,
				  patch,
				  MultiIndex<int,2>(basis1d().DeltaLmin()+1+(help/(Deltaj-2)),
						    basis1d().DeltaLmin()+1+(help%(Deltaj-2))),
				  this),
			    it->second);
	} else {
	  patch = 3+(help-3*(Deltaj-2)*(Deltaj-2))/(Deltaj-2);
	  help -= 3*(Deltaj-2)*(Deltaj-2)+(patch-3)*(Deltaj-2);

	  if (patch == 3)
	    c.add_coefficient(Index(j,
				    zero,
				    patch,
				    MultiIndex<int,2>(basis1d().DeltaLmin()+1+help,
						      0),
				    this),
			      it->second);
	  else
	    c.add_coefficient(Index(j,
				    zero,
				    patch,
				    MultiIndex<int,2>(0,
						      basis1d().DeltaLmin()+1+help),
				    this),
			      it->second);
	}
      }
    }
  }

  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::evaluate
  (const typename LDomainBasis<IntervalBasis>::Index& lambda,
   const int patch,
   const Array1D<double>& xlist,
   const Array1D<double>& ylist,
   Array1D<double>& funcvalues) const
  {
    // compute the generator expansion of psi_lambda
    InfiniteVector<double,Index> gcoeffs;
    typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()), gcoeffs_end(gcoeffs.begin());
    const int ecode = lambda.e()[0]+2*lambda.e()[1];
    const int level = lambda.j() + (ecode == 0 ? 0 : 1);
    if (ecode == 0) {
      gcoeffs.set_coefficient(lambda, 1.0);
      it = gcoeffs.begin();
      gcoeffs_end = gcoeffs.end();
    } else {
      reconstruct_1(lambda, level, it, gcoeffs_end);
    }
    
    const int nx(xlist.size());
    const int ny(ylist.size());
    Array1D<double> fx(nx);
    Array1D<double> fy(ny);

    funcvalues.resize(nx*ny);
    for (int i = 0; i < nx*ny; i++)
      funcvalues[i] = 0;

    // iterate through the generators and collect the point values on patch <patch>
    for (; it != gcoeffs_end; ++it) {
      if (it.index().p() == patch) { // i.e. lambda.p <= 2 and lambda.p == p
 	// "patch generator"
 	basis1d().evaluate(0,
 			   level, 0, it.index().k()[0],
 			   xlist,
 			   fx);
	
 	basis1d().evaluate(0,
 			   level, 0, it.index().k()[1],
 			   ylist,
 			   fy);

 	for (int n = 0, id = 0; n < nx; n++) {
  	  if (fx[n] != 0) {
  	    const double help(*it * fx[n]);
  	    for (int m = 0; m < ny; m++, id++) {
 	      funcvalues[id] += help * fy[m];
  	    }
  	  } else {
  	    id += ny;
	  }
 	}	  
      } else {
	switch(it.index().p()) {
	case 3:
	  {
	    switch(patch) {
	    case 0:
	      basis1d().evaluate(0,
				 level, 0, it.index().k()[0],
				 xlist,
				 fx);
	      
	      basis1d().evaluate(0,
				 level, 0, basis1d().DeltaLmin(),
				 ylist,
				 fy);
	    
	      for (int n = 0, id = 0; n < nx; n++) {
 		if (fx[n] != 0) {
		  const double help(M_SQRT1_2 * *it * fx[n]);
		  for (int m = 0; m < ny; m++, id++) {
		    funcvalues[id] += help * fy[m];
		  }
		} else {
		  id += ny;
		}
	      } 
	      break;
	    case 1:
	      basis1d().evaluate(0,
				 level, 0, it.index().k()[0],
				 xlist,
				 fx);
	      
	      basis1d().evaluate(0,
				 level, 0, basis1d().DeltaRmax(level),
 				 ylist,
				 fy);

	      for (int n = 0, id = 0; n < nx; n++) {
		if (fx[n] != 0) {
		  const double help(M_SQRT1_2 * *it * fx[n]);
		  for (int m = 0; m < ny; m++, id++) {
		    funcvalues[id] += help * fy[m];
		  }
		} else {
		  id += ny;
		}
	      }
	      break;
	    default:
	      break;
	    }
    	  }
	  break;
	case 4:
	  {
	    switch(patch) {
	    case 1:
	      basis1d().evaluate(0,
				 level, 0, basis1d().DeltaRmax(level),
				 xlist,
				 fx);
	      
	      basis1d().evaluate(0,
				 level, 0, it.index().k()[1],
				 ylist,
				 fy);

	      for (int n = 0, id = 0; n < nx; n++) {
		if (fx[n] != 0) {
		  const double help(M_SQRT1_2 * *it * fx[n]);
		  for (int m = 0; m < ny; m++, id++) {
		    funcvalues[id] += help * fy[m];
		  }
		} else {
		  id += ny;
		}
	      }
	      break;
	    case 2:
	      basis1d().evaluate(0,
				 level, 0, basis1d().DeltaLmin(),
				 xlist,
				 fx);
	      
	      basis1d().evaluate(0,
				 level, 0, it.index().k()[1],
				 ylist,
				 fy);
	      
	      for (int n = 0, id = 0; n < nx; n++) {
		if (fx[n] != 0) {
		  const double help(M_SQRT1_2 * *it * fx[n]);
		  for (int m = 0; m < ny; m++, id++) {
		    funcvalues[id] += help * fy[m];
		  }
		} else {
		  id += ny;
		}
	      }
	      break;
	    default:
	      break;
	    }
	  }
	  break;
	default:
	  break;
	}
      }
    }
  }

}

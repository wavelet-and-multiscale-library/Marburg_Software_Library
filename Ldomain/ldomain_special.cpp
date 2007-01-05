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

  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::reconstruct_1
  (const Index& lambda,
   const int j,
   InfiniteVector<double, Index>& c) const
  {
    assert(j == lambda.j()+1); // we only support a level difference of one
    
    typedef std::map<size_type,double> V;

    const int ecode(lambda.e()[0]+2*lambda.e()[1]);
    switch(ecode) {
    case 0:
      // generator
      switch(lambda.p()) {
      case 0: {
	// psi_lambda decomposes into generators on patch 0 alone
	// apply kron(M#,M#), cf. KroneckerMatrix::apply()
	// 1. first decide which block number and subindex corresponds to the given generator
	const size_type block_nr    = lambda.k()[0]-(basis1d().DeltaLmin()+1);
	const size_type block_index = lambda.k()[1]-(basis1d().DeltaLmin()+1);
// 	cout << "block_nr=" << block_nr << endl;
// 	cout << "block_index=" << block_index << endl;

 	// 2. get column of second factor M#
 	V z1, z2, z3;
 	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(lambda.j());
 	basis1d().Mj0_.apply_central_block(z1, z2);
//  	cout << "- column " << block_index << " of second factor M#:" << endl;
//  	for (V::const_iterator it(z2.begin()); it != z2.end(); ++it)
//  	  cout << "c[" << it->first << "]=" << it->second << endl;

 	// 3. get column of first factor M#
 	z1.clear();
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);
//  	cout << "- column " << block_nr << " of first factor M#:" << endl;
//  	for (V::const_iterator it(z3.begin()); it != z3.end(); ++it)
//  	  cout << "c[" << it->first << "]=" << it->second << endl;
	
 	// 4. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    c.add_coefficient(Index(lambda.j()+1,
				    MultiIndex<int,2>(0,0),
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
 	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(lambda.j());
 	basis1d().Mj0_.apply_central_block(z1, z2);

 	z1.clear();
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    c.add_coefficient(Index(lambda.j()+1,
				    MultiIndex<int,2>(0,0),
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
 	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(lambda.j());
 	basis1d().Mj0_.apply_central_block(z1, z2);

 	z1.clear();
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    c.add_coefficient(Index(lambda.j()+1,
				    MultiIndex<int,2>(0,0),
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
	z1[0] = 1.0;
	basis1d().Mj0_.set_level(lambda.j());
	basis1d().Mj0_.apply(z1, z2); // we have to neglect the first entry of z2 later
//   	cout << "- second factor ML:" << endl;
//   	for (V::const_iterator it(z2.begin()); it != z2.end(); ++it)
//   	  cout << "c[" << it->first << "]=" << it->second << endl;

 	// 2. get column of first factor M#
 	z1.clear();
 	const size_type block_nr = lambda.k()[0]-(basis1d().DeltaLmin()+1);
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);
//   	cout << "- column " << block_nr << " of first factor M#:" << endl;
//   	for (V::const_iterator it(z3.begin()); it != z3.end(); ++it)
//   	  cout << "c[" << it->first << "]=" << it->second << endl;

 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(++(z2.begin())); it2 != z2.end(); ++it2)
	    c.add_coefficient(Index(lambda.j()+1,
				    MultiIndex<int,2>(0,0),
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
	z1[basis1d().Deltasize(lambda.j())-1] = 1.0;
	z2.clear();
	basis1d().Mj0_.apply(z1, z2); // we have to neglect the last entry of z2 later
//    	cout << "- second factor MR:" << endl;
//    	for (V::const_iterator it(z2.begin()); it != z2.end(); ++it)
//    	  cout << "c[" << it->first << "]=" << it->second << endl;

	// 2. get column of first factor M#: this has already been done above
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::reverse_iterator it2(++(z2.rbegin())); it2 != z2.rend(); ++it2)
	    c.add_coefficient(Index(lambda.j()+1,
				    MultiIndex<int,2>(0,0),
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
				  MultiIndex<int,2>(0,0),
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
	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(lambda.j());
	basis1d().Mj0_.apply_central_block(z1, z2);
//   	cout << "- column " << block_index << " of second factor M#:" << endl;
//   	for (V::const_iterator it(z2.begin()); it != z2.end(); ++it)
//   	  cout << "c[" << it->first << "]=" << it->second << endl;

 	// 2. get column of first factor MR
	z1.clear();
	z1[basis1d().Deltasize(lambda.j())-1] = 1.0;
	basis1d().Mj0_.apply(z1, z3); // we have to neglect the last entry of z3 later
	
	// 3. combine results
 	for (V::reverse_iterator it3(++(z3.rbegin())); it3 != z3.rend(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    c.add_coefficient(Index(lambda.j()+1,
				    MultiIndex<int,2>(0,0),
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
	z1[0] = 1.0;
	z3.clear();
	basis1d().Mj0_.apply(z1, z3); // we have to neglect the first entry of z3 later

	// 3. combine results
 	for (V::const_iterator it3(++(z3.begin())); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    c.add_coefficient(Index(lambda.j()+1,
				    MultiIndex<int,2>(0,0),
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
				  MultiIndex<int,2>(0,0),
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
      break;
    case 2:
      // (0,1)-wavelet
      break;
    case 3:
      // (1,1)-wavelet
      break;
    default:
      break;
    }
  }


}

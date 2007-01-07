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

 	// 2. get column of second factor M#
 	V z1, z2, z3;
 	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(lambda.j());
 	basis1d().Mj0_.apply_central_block(z1, z2);

 	// 3. get column of first factor M#
 	z1.clear();
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);
	
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

 	// 2. get column of first factor M#
 	z1.clear();
 	const size_type block_nr = lambda.k()[0]-(basis1d().DeltaLmin()+1);
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);

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
    for (V::const_iterator it(x.begin()); it != x.end(); ++it) {
      size_type help(it->first);
      unsigned int patch = 0;
      if (help < 3*(Deltaj-2)*(Deltaj-2)) {
	patch = help / ((Deltaj-2)*(Deltaj-2));
	help -= patch * (Deltaj-2)*(Deltaj-2);
      } else {
	patch = 3+(help-3*(Deltaj-2)*(Deltaj-2))/(Deltaj-2);
	help -= 3*(Deltaj-2)*(Deltaj-2)+(patch-3)*(Deltaj-2);
      }
      
      switch(patch) {
      case 0:
      case 1:
      case 2:
	c.add_coefficient(Index(j,
				MultiIndex<int,2>(0,0),
				patch,
				MultiIndex<int,2>(basis1d().DeltaLmin()+1+(help/(Deltaj-2)),
						  basis1d().DeltaLmin()+1+(help%(Deltaj-2))),
				this),
			  it->second);
	break;
      case 3:
	c.add_coefficient(Index(j,
				MultiIndex<int,2>(0,0),
				patch,
				MultiIndex<int,2>(basis1d().DeltaLmin()+1+help,
						  0),
				this),
			  it->second);
	break;
      case 4:
	c.add_coefficient(Index(j,
				MultiIndex<int,2>(0,0),
				patch,
				MultiIndex<int,2>(0,
						  basis1d().DeltaLmin()+1+help),
				this),
			  it->second);
	break;
      default:
	break;
      }
    }
  }

  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::apply_Mj0
  (const int j,
   const std::map<size_type,double>& x, 
   std::map<size_type,double>& y) const
  {
    typedef std::map<size_type,double> V;
    
    const unsigned int Deltaj   = basis1d().Deltasize(j);
    const unsigned int Deltajp1 = basis1d().Deltasize(j+1);
    
    for (V::const_iterator itx(x.begin()); itx != x.end(); ++itx) {
      // determine patch number
      const unsigned int patch =
	itx->first < 3*(Deltaj-2)*(Deltaj-2)
	? itx->first / ((Deltaj-2)*(Deltaj-2))
	: 3+(itx->first-3*(Deltaj-2)*(Deltaj-2))/(Deltaj-2);
      
      switch(patch) {
      case 0: {
	// psi_lambda decomposes into generators on patch 0 alone
	// apply kron(M#,M#), cf. KroneckerMatrix::apply()
 	// 1. get column of second factor M#
 	V z1, z2, z3;
	const size_type block_index = itx->first % (Deltaj-2);
 	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(j);
 	basis1d().Mj0_.apply_central_block(z1, z2);
	
 	// 2. get column of first factor M#
 	z1.clear();
	const size_type block_nr    = itx->first / (Deltaj-2);
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);
	
 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
	break;
      case 1: {
	// psi_lambda decomposes into generators on patch 1 alone
	// apply kron(M#,M#)
	const size_type block_nr    = (itx->first-(Deltaj-2)*(Deltaj-2)) / (Deltaj-2);
	const size_type block_index = (itx->first-(Deltaj-2)*(Deltaj-2)) % (Deltaj-2);
	
 	V z1, z2, z3;
 	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(j);
 	basis1d().Mj0_.apply_central_block(z1, z2);

 	z1.clear();
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
      }
      case 2: {
	// psi_lambda decomposes into generators on patch 2 alone
	// apply kron(M#,M#)
	const size_type block_nr    = (itx->first-2*(Deltaj-2)*(Deltaj-2)) / (Deltaj-2);
	const size_type block_index = (itx->first-2*(Deltaj-2)*(Deltaj-2)) % (Deltaj-2);

 	V z1, z2, z3;
 	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(j);
 	basis1d().Mj0_.apply_central_block(z1, z2);

 	z1.clear();
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first]
	      += itx->second * it2->second * it3->second;
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
	basis1d().Mj0_.set_level(j);
	basis1d().Mj0_.apply(z1, z2); // we have to neglect the first entry of z2 later

 	// 2. get column of first factor M#
 	z1.clear();
 	const size_type block_nr = itx->first-3*(Deltaj-2)*(Deltaj-2);
 	z1[block_nr] = 1.0;
 	basis1d().Mj0_.apply_central_block(z1, z3);

 	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(++(z2.begin())); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltajp1-2)+it2->first-1] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

	// contribution from patch 1
	//
	// apply kron(M#,MR):
 	// 1. get column of second factor MR
	z1.clear();
	z1[basis1d().Deltasize(j)-1] = 1.0;
	z2.clear();
	basis1d().Mj0_.apply(z1, z2); // we have to neglect the last entry of z2 later

	// 2. get column of first factor M#: this has already been done above
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::reverse_iterator it2(++(z2.rbegin())); it2 != z2.rend(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+it3->first*(Deltajp1-2)+it2->first-1] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

	// contribution from patch 3
	//
	// apply kron(M#,Mtopleft):
	// 1. get column of first factor M#: this has already been done above

	// 2. get top left entry of Mj0
	const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	
	// 3. combine results
	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	  y[3*(Deltajp1-2)*(Deltajp1-2)+it3->first]
	    += itx->second * it3->second * Mtopleft;
      }
	break;
      case 4: {
	// psi_lambda decomposes into generators on patches 1,2 and 4

	// contribution from patch 1
	//
	// apply kron(MR,M#):
 	// 1. get column of second factor M#
	V z1, z2, z3;
	const size_type block_index = itx->first-(3*(Deltaj-2)+1)*(Deltaj-2);
	z1[block_index] = 1.0;
	basis1d().Mj0_.set_level(j);
	basis1d().Mj0_.apply_central_block(z1, z2);

 	// 2. get column of first factor MR
	z1.clear();
	z1[basis1d().Deltasize(j)-1] = 1.0;
	basis1d().Mj0_.apply(z1, z3); // we have to neglect the last entry of z3 later
	
	// 3. combine results
 	for (V::reverse_iterator it3(++(z3.rbegin())); it3 != z3.rend(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltajp1-2)*(Deltajp1-2)+(it3->first-1)*(Deltajp1-2)+it2->first] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

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
	    y[2*(Deltajp1-2)*(Deltajp1-2)+(it3->first-1)*(Deltajp1-2)+it2->first] // note the "-1"
	      += itx->second * it2->second * it3->second * M_SQRT1_2;

	// contribution from patch 4
	//
	// apply kron(Mtopleft,M#):
	// 1. get column of second factor M#: this has already been done above

	// 2. get top left entry of Mj0
	const double Mtopleft = basis1d().Mj0_.get_entry(0,0);
	
	// 3. combine results
	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	  y[(3*(Deltajp1-2)+1)*(Deltajp1-2)+it2->first]
	    += itx->second * it2->second * Mtopleft;
      }
	break;
      default:
	break;
      }
    }
  }
  
  template <int d, int dT>
  void
  LDomainBasis<SplineBasis<d,dT,DS_construction> >::apply_Mj0T_transposed
  (const int j,
   const std::map<size_type,double>& x, 
   std::map<size_type,double>& y) const
  {
    typedef std::map<size_type,double> V;
    
    const unsigned int Deltaj   = basis1d().Deltasize(j);
    const unsigned int Deltajp1 = basis1d().Deltasize(j+1);
    
    for (V::const_iterator itx(x.begin()); itx != x.end(); ++itx) {
//       cout << "apply_Mj0T_t(): number " << itx->first << ",  value " << itx->second << endl;

      // determine patch number
      const unsigned int patch =
	itx->first < 3*(Deltajp1-2)*(Deltajp1-2)
	? itx->first / ((Deltajp1-2)*(Deltajp1-2))
	: 3+(itx->first-3*(Deltajp1-2)*(Deltajp1-2))/(Deltajp1-2);
//       cout << "apply_Mj0T_t(): patch number " << patch << endl;
      
      switch(patch) {
      case 0: {
	// apply kron((M#)^T,(M#)^T):
 	// 1. get column of second factor (M#)^T
 	V z1, z2, z3;
	const size_type block_index = itx->first % (Deltajp1-2);
 	z1[block_index] = 1.0;
	basis1d().Mj0T_.set_level(j);
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);
// 	cout << "row " << block_index << " of Mj0T is" << endl;
//  	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
//  	  cout << it2->first << ": " << it2->second << endl;

 	// 2. get column of first factor (M#)^T
 	z1.clear();
	const size_type block_nr    = itx->first / (Deltajp1-2);
 	z1[block_nr] = 1.0;
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z3);
// 	cout << "row " << block_nr << " of Mj0T is" << endl;
//  	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
//  	  cout << it3->first << ": " << it3->second << endl;
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[it3->first*(Deltaj-2)+it2->first]
	      += itx->second * it2->second * it3->second;

 	// apply kron((M#)^T,(ML)^T)
	// 1. get column of second factor (ML)^T (which is a number)
	const double ML_entry = basis1d().Mj0T_.get_entry(block_index+1, 0);

	// 2. get column of first factor (M#)^T: this has already been done above

	// 3. combine results
	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	  y[3*(Deltaj-2)*(Deltaj-2)+it3->first]
	    += itx->second * ML_entry * it3->second * M_SQRT1_2;
      }
	break;
      case 1: {
	// apply kron((M#)^T,(M#)^T):
 	// 1. get column of second factor (M#)^T
 	V z1, z2, z3;
	const size_type block_index = (itx->first-(Deltajp1-2)*(Deltajp1-2)) % (Deltajp1-2);
 	z1[block_index] = 1.0;
	basis1d().Mj0T_.set_level(j);
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

 	// 2. get column of first factor (M#)^T
 	z1.clear();
	const size_type block_nr    = (itx->first-(Deltajp1-2)*(Deltajp1-2)) / (Deltajp1-2);
 	z1[block_nr] = 1.0;
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z3);
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[(Deltaj-2)*(Deltaj-2)+it3->first*(Deltaj-2)+it2->first]
	      += itx->second * it2->second * it3->second;

	// apply kron((M#)^T,(MR)^T)
	// 1. get column of second factor (MR)^T (which is a number)
	const double MR_entry = basis1d().Mj0T_.get_entry(block_index+1, Deltaj-1);

	// 2. get column of first factor (M#)^T: this has already been done above
	
	// 3. combine results
	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
	  y[3*(Deltaj-2)*(Deltaj-2)+it3->first]
	    += itx->second * MR_entry * it3->second * M_SQRT1_2;

	// apply kron((MR)^T,(M#)^T)
	// 1. get column of second factor (M#)^T: this has already been done above

	// 2. get column of first factor (MR)^T (which is a number)
	const double MR_entry2 = basis1d().Mj0T_.get_entry(block_nr+1, Deltaj-1);
	
	// 3. combine results
	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	  y[(3*(Deltaj-2)+1)*(Deltaj-2)+it2->first]
	    += itx->second * MR_entry2 * it2->second * M_SQRT1_2;
      }
	break;
      case 2: {
	// apply kron((M#)^T,(M#)^T):
 	// 1. get column of second factor (M#)^T
 	V z1, z2, z3;
	const size_type block_index = (itx->first-2*(Deltajp1-2)*(Deltajp1-2)) % (Deltajp1-2);
 	z1[block_index] = 1.0;
	basis1d().Mj0T_.set_level(j);
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

 	// 2. get column of first factor (M#)
 	z1.clear();
	const size_type block_nr    = (itx->first-2*(Deltajp1-2)*(Deltajp1-2)) / (Deltajp1-2);
 	z1[block_nr] = 1.0;
 	basis1d().Mj0T_.apply_central_block_transposed(z1, z3);
	
	// 3. combine results
 	for (V::const_iterator it3(z3.begin()); it3 != z3.end(); ++it3)
 	  for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	    y[2*(Deltaj-2)*(Deltaj-2)+it3->first*(Deltaj-2)+it2->first]
	      += itx->second * it2->second * it3->second;

	// apply kron((ML)^T,(M#)^T)
	// 1. get column of second factor (M#)^T: this has already been done above
	
	// 2. get column of first factor (ML)^T (which is a number)
	const double ML_entry = basis1d().Mj0T_.get_entry(block_nr+1, 0);

	// 3. combine results
	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
	  y[(3*(Deltaj-2)+1)*(Deltaj-2)+it2->first]
	    += itx->second * ML_entry * it2->second * M_SQRT1_2;	
      }
	break;
      case 3: {
	// apply kron(Mtopleft,(M#)^T):
	// 1. get column of second factor (M#)^T
	V z1, z2;
 	const size_type block_nr = itx->first-3*(Deltajp1-2)*(Deltajp1-2);
	z1[block_nr] = 1.0;
	basis1d().Mj0T_.set_level(j);
	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

	// 2. get top left entry of Mj0T
	const double Mtopleft = basis1d().Mj0T_.get_entry(0,0);

	// 3. combine the results
 	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
 	  y[3*(Deltaj-2)*(Deltaj-2)+it2->first]
 	    += itx->second * it2->second * Mtopleft;
      }
	break;
      case 4: {
	// apply kron(Mtopleft,(M#)^T):
	V z1, z2;
 	const size_type block_nr = itx->first-(3*(Deltajp1-2)+1)*(Deltajp1-2);
	z1[block_nr] = 1.0;
	basis1d().Mj0T_.set_level(j);
	basis1d().Mj0T_.apply_central_block_transposed(z1, z2);

	// 2. get top left entry of Mj0T
	const double Mtopleft = basis1d().Mj0T_.get_entry(0,0);

	// 3. combine the results
 	for (V::const_iterator it2(z2.begin()); it2 != z2.end(); ++it2)
 	  y[(3*(Deltaj-2)+1)*(Deltaj-2)+it2->first]
 	    += itx->second * it2->second * Mtopleft;	
      }
	break;
      default:
	break;
      }
    }
  }
}

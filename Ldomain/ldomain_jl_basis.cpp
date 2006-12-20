// implementation for ldomain_jl_basis.h

#include <map>
#include <utils/fixed_array1d.h>

using MathTL::FixedArray1D;

namespace WaveletTL
{
  LDomainJLBasis::LDomainJLBasis()
    : j0_(1)
  {
  }
  
  void
  LDomainJLBasis::reconstruct_1_1d(const RMWIndex& lambda, const int j,
				   InfiniteVector<double,RMWIndex>& c) const {
    if (lambda.j() >= j) {
      c.add_coefficient(lambda, 1.0);
    } else {
      assert(lambda.j()+1 == j); // for the moment we only allow one step at a time
      
      if (lambda.e() == 0) {
	// generator
	if (lambda.c() == 0) {
	  // type phi_0
	  // phi_0(x) = 1/2*phi_0(2*x+1)+phi_0(2*x)+1/2*phi_0(2*x-1)+3/4*phi_1(2*x+1)-3/4*phi_1(2*x-1)
	  
	  int m = 2*lambda.k()-1; // m-2k=-1
	  c.add_coefficient(RMWIndex(j, 0, 0, m), 0.5*M_SQRT1_2);  // phi_0(2x+1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), 0.75*M_SQRT1_2); // phi_1(2x+1)
	  
	  // m=2k <-> m-2k=0
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 0, m), M_SQRT1_2); // phi_0(2x)
	  
	  // m=2k+1 <-> m-2k=1
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 0, m), 0.5*M_SQRT1_2);   // phi_0(2x-1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), -0.75*M_SQRT1_2); // phi_1(2x-1)
	} else {
	  // lambda.c() == 1
	  // phi_1(x) = -1/8*phi_0(2*x+1)+1/8*phi_0(2*x-1)-1/8*phi_1(2*x+1)+1/2*phi_1(2*x)-1/8*phi_1(2*x-1)
	  
	  int m = 2*lambda.k()-1; // m-2k=-1
	  c.add_coefficient(RMWIndex(j, 0, 0, m), -0.125*M_SQRT1_2); // phi_0(2x+1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), -0.125*M_SQRT1_2); // phi_1(2x+1)
	  
	  // m=2k <-> m-2k=0
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 1, m), 0.5*M_SQRT1_2); // phi_1(2x)
	  
	  // m=2k+1 <-> m-2k=1
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 0, m), 0.125*M_SQRT1_2);  // phi_0(2x-1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), -0.125*M_SQRT1_2); // phi_1(2x-1)
	}
      } else { // lambda.e() == 1
	// wavelet
	if (lambda.c() == 0) {
	  // type psi_0
	  // psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)
	  
	  int m = 2*lambda.k()-1; // m-2k=-1
	  c.add_coefficient(RMWIndex(j, 0, 0, m), -2.0*M_SQRT1_2);  // phi_0(2x+1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), -21.0*M_SQRT1_2); // phi_1(2x+1)
	  
	  // m=2k <-> m-2k=0
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 0, m), 4.0*M_SQRT1_2); // phi_0(2x)
	  
	  // m=2k+1 <-> m-2k=1
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 0, m), -2.0*M_SQRT1_2); // phi_0(2x-1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), 21.0*M_SQRT1_2); // phi_1(2x-1)
	} else {
	  // type psi_1
	  // psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)
	  
	  int m = 2*lambda.k()-1; // m-2k=-1
	  c.add_coefficient(RMWIndex(j, 0, 0, m), M_SQRT1_2);     // phi_0(2x+1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), 9.0*M_SQRT1_2); // phi_1(2x+1)
	  
	  // m=2k <-> m-2k=0
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 1, m), 12.0*M_SQRT1_2); // phi_1(2x)
	  
	  // m=2k+1 <-> m-2k=1
	  m++;
	  c.add_coefficient(RMWIndex(j, 0, 0, m), -M_SQRT1_2);    // phi_0(2x-1)
	  c.add_coefficient(RMWIndex(j, 0, 1, m), 9.0*M_SQRT1_2); // phi_1(2x-1)
	}
      }
    }
  }

  void
  LDomainJLBasis::reconstruct_1(const Index& lambda,
				const int j,
				InfiniteVector<double, Index>& c) const {
    if (lambda.j() >= j) {
      c.add_coefficient(lambda, 1.0);
    } else {
      assert(lambda.j()+1 == j); // we assume that we only have to perform one step
      
#if 0
      // first variant:
      // perform 2 1D reconstruct_1() calls, only store the k parameters
      FixedArray1D<InfiniteVector<double,int>,2> coeffs0, coeffs1; // for c=0,1 resp.
      for (int n = 0; n <= 1; n++) {
	if (lambda.e()[n] == 0) {
	  // generator
	  if (lambda.c()[n] == 0) {
	    // type phi_0
	    // phi_0(x) = 1/2*phi_0(2*x+1)+phi_0(2*x)+1/2*phi_0(2*x-1)+3/4*phi_1(2*x+1)-3/4*phi_1(2*x-1)
	    
	    int m = 2*lambda.k()[n]-1; // m-2k=-1
	    coeffs0[n].add_coefficient(m, 0.5*M_SQRT1_2);  // phi_0(2x+1)
	    coeffs1[n].add_coefficient(m, 0.75*M_SQRT1_2); // phi_1(2x+1)
	    
	    // m=2k <-> m-2k=0
	    m++;
	    coeffs0[n].add_coefficient(m, M_SQRT1_2); // phi_0(2x)
	    
	    // m=2k+1 <-> m-2k=1
	    m++;
	    coeffs0[n].add_coefficient(m, 0.5*M_SQRT1_2);   // phi_0(2x-1)
	    coeffs1[n].add_coefficient(m, -0.75*M_SQRT1_2); // phi_1(2x-1)
	  } else {
	    // lambda.c()[n] == 1
	    // phi_1(x) = -1/8*phi_0(2*x+1)+1/8*phi_0(2*x-1)-1/8*phi_1(2*x+1)+1/2*phi_1(2*x)-1/8*phi_1(2*x-1)

	    int m = 2*lambda.k()[n]-1; // m-2k=-1
	    coeffs0[n].add_coefficient(m, -0.125*M_SQRT1_2); // phi_0(2x+1)
	    coeffs1[n].add_coefficient(m, -0.125*M_SQRT1_2); // phi_1(2x+1)

	    // m=2k <-> m-2k=0
	    m++;
	    coeffs1[n].add_coefficient(m, 0.5*M_SQRT1_2); // phi_1(2x)

	    // m=2k+1 <-> m-2k=1
	    m++;
	    coeffs0[n].add_coefficient(m, 0.125*M_SQRT1_2);  // phi_0(2x-1)
	    coeffs1[n].add_coefficient(m, -0.125*M_SQRT1_2); // phi_1(2x-1)
	  }
	} else { // lambda.e() == 1
	  // wavelet
	  if (lambda.c()[n] == 0) {
	    // type psi_0
	    // psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)

	    int m = 2*lambda.k()[n]-1; // m-2k=-1
	    coeffs0[n].add_coefficient(m, -2.0*M_SQRT1_2);  // phi_0(2x+1)
	    coeffs1[n].add_coefficient(m, -21.0*M_SQRT1_2); // phi_1(2x+1)

	    // m=2k <-> m-2k=0
	    m++;
	    coeffs0[n].add_coefficient(m, 4.0*M_SQRT1_2); // phi_0(2x)

	    // m=2k+1 <-> m-2k=1
	    m++;
	    coeffs0[n].add_coefficient(m, -2.0*M_SQRT1_2); // phi_0(2x-1)
	    coeffs1[n].add_coefficient(m, 21.0*M_SQRT1_2); // phi_1(2x-1)
	  } else {
	    // type psi_1
	    // psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)

	    int m = 2*lambda.k()[n]-1; // m-2k=-1
	    coeffs0[n].add_coefficient(m, M_SQRT1_2);     // phi_0(2x+1)
	    coeffs1[n].add_coefficient(m, 9.0*M_SQRT1_2); // phi_1(2x+1)

	    // m=2k <-> m-2k=0
	    m++;
	    coeffs1[n].add_coefficient(m, 12.0*M_SQRT1_2); // phi_1(2x)

	    // m=2k+1 <-> m-2k=1
	    m++;
	    coeffs0[n].add_coefficient(m, -M_SQRT1_2);    // phi_0(2x-1)
	    coeffs1[n].add_coefficient(m, 9.0*M_SQRT1_2); // phi_1(2x-1)
	  }
	}
      } // end for (int n = 0; n <= 1; n++)
      
      // directly add all tensor product generators involved
      for (InfiniteVector<double,int>::const_iterator it00 = coeffs0[0].begin();
	   it00 != coeffs0[0].end(); ++it00) {
	for (InfiniteVector<double,int>::const_iterator it01 = coeffs0[1].begin();
	     it01 != coeffs0[1].end(); ++it01) {
	  // c=(0,0)
	  if (index_is_valid(lambda.j()+1, 0, 0, 0, 0, it00.index(), it01.index()))
	    c.add_coefficient(Index(lambda.j()+1,
				    Index::type_type(0,0),
				    Index::component_type(0,0),
				    Index::translation_type(it00.index(), it01.index())),
			      *it00 * *it01);  
	}
	for (InfiniteVector<double,int>::const_iterator it11 = coeffs1[1].begin();
	     it11 != coeffs1[1].end(); ++it11) {
	  // c=(0,1)
	  if (index_is_valid(lambda.j()+1, 0, 0, 0, 1, it00.index(), it11.index()))
	    c.add_coefficient(Index(lambda.j()+1,
				    Index::type_type(0,0),
				    Index::component_type(0,1),
				    Index::translation_type(it00.index(), it11.index())),
			      *it00 * *it11);  
	}
      }
      for (InfiniteVector<double,int>::const_iterator it10 = coeffs1[0].begin();
	   it10 != coeffs1[0].end(); ++it10) {
	for (InfiniteVector<double,int>::const_iterator it01 = coeffs0[1].begin();
	     it01 != coeffs0[1].end(); ++it01) {
	  // c=(1,0)
	  if (index_is_valid(lambda.j()+1, 0, 0, 1, 0, it10.index(), it01.index()))
	    c.add_coefficient(Index(lambda.j()+1,
				    Index::type_type(0,0),
				    Index::component_type(1,0),
				    Index::translation_type(it10.index(), it01.index())),
			      *it10 * *it01);  
	}
	for (InfiniteVector<double,int>::const_iterator it11 = coeffs1[1].begin();
	     it11 != coeffs1[1].end(); ++it11) {
	  // c=(1,1)
	  if (index_is_valid(lambda.j()+1, 0, 0, 1, 1, it10.index(), it11.index()))
	    c.add_coefficient(Index(lambda.j()+1,
				    Index::type_type(0,0),
				    Index::component_type(1,1),
				    Index::translation_type(it10.index(), it11.index())),
			      *it10 * *it11);  
	}
      }
#else
      // second variant:
      // perform 2 1D reconstruct_1_1d() calls
      InfiniteVector<double,RMWIndex> coeffs0, coeffs1;
      reconstruct_1_1d(RMWIndex(lambda.j(),lambda.e()[0],lambda.c()[0],lambda.k()[0]),lambda.j()+1,coeffs0);
      reconstruct_1_1d(RMWIndex(lambda.j(),lambda.e()[1],lambda.c()[1],lambda.k()[1]),lambda.j()+1,coeffs1);

      // directly add all tensor product generators involved
      for (InfiniteVector<double,RMWIndex>::const_iterator it0 = coeffs0.begin();
	   it0 != coeffs0.end(); ++it0) {
	for (InfiniteVector<double,RMWIndex>::const_iterator it1 = coeffs1.begin();
	     it1 != coeffs1.end(); ++it1) {
	  if (index_is_valid(lambda.j()+1, 0, 0, it0.index().c(), it1.index().c(), it0.index().k(), it1.index().k()))
	    c.add_coefficient(Index(lambda.j()+1,
				    Index::type_type(0,0),
				    Index::component_type(it0.index().c(), it1.index().c()),
				    Index::translation_type(it0.index().k(), it1.index().k())),
			      *it0 * *it1);  
	}
      }
#endif
    }
  }

  LDomainJLIndex
  LDomainJLBasis::first_generator(const int j) const
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e;      // e=(0,0)
    const LDomainJLIndex::component_type c; // c=(0,0)
    const LDomainJLIndex::translation_type k(1-(1<<j), 1-(1<<j));
    return LDomainJLIndex(j, e, c, k);
  }

  LDomainJLIndex
  LDomainJLBasis::last_generator(const int j) const
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e; // e=(0,0)
    const LDomainJLIndex::component_type c(1,1);
    const LDomainJLIndex::translation_type k(1<<j,0); // generator sits at upper right corner of the "forefoot"
    return LDomainJLIndex(j, e, c, k);
  }
  
  LDomainJLIndex
  LDomainJLBasis::first_wavelet(const int j) const
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e(0,1);
    const LDomainJLIndex::component_type c(0,0);
    const LDomainJLIndex::translation_type k(1-(1<<j),1-(1<<j));
    return LDomainJLIndex(j, e, c, k);
  }

  LDomainJLIndex
  LDomainJLBasis::first_wavelet(const int j,
				const LDomainJLIndex::type_type& e) const
  {
    assert(j >= 1);
    const LDomainJLIndex::component_type c(0,0);
    const LDomainJLIndex::translation_type k(1-(1<<j),1-(1<<j));
    return LDomainJLIndex(j, e, c, k);
  }
  

  LDomainJLIndex
  LDomainJLBasis::last_wavelet(const int j) const
  {
    assert(j >= 1);
    const LDomainJLIndex::type_type e(1,1);
    const LDomainJLIndex::component_type c(1,1);
    const LDomainJLIndex::translation_type k(1<<j,0); // wavelet sits at upper right corner of the "forefoot"
    return LDomainJLIndex(j, e, c, k);
  }

  
}

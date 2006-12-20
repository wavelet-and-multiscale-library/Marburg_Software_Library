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
  LDomainJLBasis::reconstruct_1(const Index& lambda,
				const int j,
				InfiniteVector<double, Index>& c) const {
    assert(lambda.j()+1 == j); // we assume that we only have to perform one step
    
    // perform 2 1D reconstruct_1() calls, we only store the k parameters
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

//       cout << "In reconstruct_1() for lambda=" << lambda << "," << endl
// 	   << "coeffs0[0]=" << endl << coeffs0[0] << endl
// 	   << "coeffs1[0]=" << endl << coeffs1[0] << endl
// 	   << "coeffs0[1]=" << endl << coeffs0[1] << endl
// 	   << "coeffs1[1]=" << endl << coeffs1[1] << endl;
      
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
  }
  
}

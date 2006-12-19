// implementation for ldomain_jl_basis.h

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
    if (lambda.j() >= j) {
      c.add_coefficient(lambda, 1.0);
    } else {
      // perform 2 1D reconstruct_1() calls, we only store the k parameters
      FixedArray1D<InfiniteVector<double,int>,2> coeffs;
      for (int n = 0; n <= 1; n++) {
	if (lambda.e()[n] == 0) {
	  // generator
	  if (lambda.c()[n] == 0) {
	    // type phi_0
	    // we know that in the n-th direction, psi_lambda does not intersect a boundary

	    // hier geht es weiter!!!
	  }
	} 
      }





//       if (lambda.e() == 0) {
//  	// generator
// 	if (lambda.c() == 0) {
//  	  // type phi_0
//  	  // phi_0(x) = 1/2*phi_0(2*x+1)+phi_0(2*x)+1/2*phi_0(2*x-1)+3/4*phi_1(2*x+1)-3/4*phi_1(2*x-1)

// 	  int m = 2*lambda.k()-1; // m-2k=-1
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(0.5*M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(0.75*M_SQRT1_2, dhelp);
// 	  }

// 	  // m=2k <-> m-2k=0
// 	  m++;
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(M_SQRT1_2, dhelp);
// 	  }
	  
// 	  // m=2k+1 <-> m-2k=1
// 	  m++;
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))) { // phi_0(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(0.5*M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))-1) { // phi_1(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(-0.75*M_SQRT1_2, dhelp);
// 	  }
// 	} else { // lambda.c() == 1
//   	  // type phi_1
//  	  // phi_1(x) = -1/8*phi_0(2*x+1)+1/8*phi_0(2*x-1)-1/8*phi_1(2*x+1)+1/2*phi_1(2*x)-1/8*phi_1(2*x-1)
	  
// 	  int m = 2*lambda.k()-1; // m-2k=-1
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(-0.125*M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(-0.125*M_SQRT1_2, dhelp);
// 	  }

// 	  // m=2k <-> m-2k=0
// 	  m++;
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(0.5*M_SQRT1_2, dhelp);
// 	  }

// 	  // m=2k+1 <-> m-2k=1
// 	  m++;
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(0.125*M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(-0.125*M_SQRT1_2, dhelp);
// 	  }
// 	} // end if (lambda.c() == 0)
//       } else { // lambda.e() == 1
// 	// wavelet
// 	if (lambda.c() == 0) {
//  	  // type psi_0
//  	  // psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)

// 	  int m = 2*lambda.k()-1; // m-2k=-1
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(-2.0*M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(-21.0*M_SQRT1_2, dhelp);
// 	  }

// 	  // m=2k <-> m-2k=0
// 	  m++;
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(4.0*M_SQRT1_2, dhelp);
// 	  }

// 	  // m=2k+1 <-> m-2k=1
// 	  m++;
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(-2.0*M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(21.0*M_SQRT1_2, dhelp);
// 	  }
// 	} else { // lambda.c() == 1
//  	  // type psi_1
//  	  // psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)

// 	  int m = 2*lambda.k()-1; // m-2k=-1
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x+1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(9.0*M_SQRT1_2, dhelp);
// 	  }

// 	  // m=2k <-> m-2k=0
// 	  m++;
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(12.0*M_SQRT1_2, dhelp);
// 	  }

// 	  // m=2k+1 <-> m-2k=1
// 	  m++;
// 	  if (m >= 1 && m <= (1<<(lambda.j()+1))-1) { // phi_0(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 0, m), j, dhelp);
// 	    c.add(-M_SQRT1_2, dhelp);
// 	  }
// 	  if (m >= 0 && m <= (1<<(lambda.j()+1))) { // phi_1(2x-1)
// 	    InfiniteVector<double, Index> dhelp;
// 	    reconstruct_1(Index(lambda.j()+1, 0, 1, m), j, dhelp);
// 	    c.add(9.0*M_SQRT1_2, dhelp);
// 	  }
// 	}



















      
//       cout << "In reconstruct_1(), coeffs[0]=" << endl << coeffs[0] << endl
// 	   << "coeffs[1]=" << endl << coeffs[1] << endl;
    }
  }


//   template <class IBASIS, unsigned int DIM>
//   void
//   CubeBasis<IBASIS,DIM>::reconstruct_1(const Index& lambda,
// 				       const int j,
// 				       InfiniteVector<double, Index>& c) const {
//     if (lambda.j() >= j) {
//       // then we can just copy \psi_\lambda
//       c.add_coefficient(lambda, 1.0);
//     } else {
//       // reconstruct by recursion

//       typedef typename IBASIS::Index IIndex;
//       FixedArray1D<InfiniteVector<double,IIndex>,DIM> coeffs;
//       for (unsigned int i = 0; i < DIM; i++)
// 	bases_[i]->reconstruct_1(IIndex(lambda.j(), lambda.e()[i], lambda.k()[i], bases_[i]), lambda.j()+1, coeffs[i]);
      
//       if (DIM == 2) { // TODO: template specialization here!!!
// 	// directly add all tensor product wavelets needed
//  	if (lambda.j()+1 >= j) {
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
// 		 it2 != it2end; ++it2) {
// 	      c.add_coefficient(Index(lambda.j()+1,
// 				      typename Index::type_type(it1.index().e(), it2.index().e()),
// 				      typename Index::translation_type(it1.index().k(), it2.index().k()),
// 				      this),
// 				*it1 * *it2);
// 	    }
//  	} else {
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(coeffs[0].begin()), it1end(coeffs[0].end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(coeffs[1].begin()), it2end(coeffs[1].end());
// 		 it2 != it2end; ++it2) {
// 	      InfiniteVector<double,Index> d;
// 	      reconstruct_1(Index(lambda.j()+1,
// 				  typename Index::type_type(it1.index().e(), it2.index().e()),
// 				  typename Index::translation_type(it1.index().k(), it2.index().k()),
// 				  this),
// 			    j, d);
// 	      c.add(*it1 * *it2, d);
// 	    }
//  	}
//       } else {
// 	// prepare all tensor product wavelet indices needed (+values)
// 	typedef std::list<std::pair<FixedArray1D<IIndex,DIM>,double> > list_type;
// 	list_type indices;
// 	FixedArray1D<IIndex,DIM> helpindex;
// 	for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[0].begin()), itvend(coeffs[0].end());
// 	     itv != itvend; ++itv) {
// 	  helpindex[0] = itv.index();
// 	  indices.push_back(std::make_pair(helpindex, *itv));
// 	}
// 	for (unsigned int i(1); i < DIM; i++) {
// 	  list_type sofar;
// 	  sofar.swap(indices);
// 	  for (typename list_type::const_iterator it(sofar.begin()), itend(sofar.end());
// 	       it != itend; ++it) {
// 	    helpindex = it->first;
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator itv(coeffs[i].begin()), itvend(coeffs[i].end());
// 		 itv != itvend; ++itv) {
// 	      helpindex[i] = itv.index();
// 	      indices.push_back(std::make_pair(helpindex, *itv * it->second));
// 	    }
// 	  }
// 	}
	
// 	// write results into c
// 	typename Index::type_type help_e;
// 	typename Index::translation_type help_k;
//  	if (lambda.j()+1 >= j) {
// 	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
// 	       it != itend; ++it) {
// 	    for (unsigned int i = 0; i < DIM; i++) {
// 	      help_e[i] = it->first[i].e();
// 	      help_k[i] = it->first[i].k();
// 	    }
// 	    c.add_coefficient(Index(lambda.j()+1, help_e, help_k, this),
// 			      it->second);
// 	  }
// 	} else {
// 	  for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
// 	       it != itend; ++it) {
// 	    InfiniteVector<double,Index> d;
// 	    for (unsigned int i = 0; i < DIM; i++) {
// 	      help_e[i] = it->first[i].e();
// 	      help_k[i] = it->first[i].k();
// 	    }
// 	    reconstruct_1(Index(lambda.j()+1, help_e, help_k, this), j, d);
// 	    c.add(it->second, d);
// 	  }
// 	}
//       }
//     }
//   }
  

}

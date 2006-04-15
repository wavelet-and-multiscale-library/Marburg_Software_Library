// implementation for ldomain_basis.h

#include <cmath>

namespace WaveletTL
{
  template <class IBASIS>
  LDomainBasis<IBASIS>::LDomainBasis()
    : basis1d_(false, false)
  {
  }

  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct(const InfiniteVector<double, Index>& c,
				    const int j,
				    InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_1(it.index(), j, help);
      d.add(*it, help);
    }
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_t(const InfiniteVector<double, Index>& c,
				      const int j,
				      InfiniteVector<double, Index>& d) const {
    for (typename InfiniteVector<double, Index>::const_iterator it(c.begin()), itend(c.end());
	 it != itend; ++it) {
      InfiniteVector<double, Index> help;
      reconstruct_t_1(it.index(), j, help);
      d.add(*it, help);
    }
  }

  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_1(const Index& lambda,
				      const int j,
				      InfiniteVector<double, Index>& c) const {
    typedef typename IBASIS::Index IIndex;
    
    if (lambda.j() >= j) {
      // then we can just copy psi_lambda
      c.add_coefficient(lambda, 1.0);
    } else {
      // For the reconstruction of psi_lambda, we have to compute
      // the corresponding column of the transformation matrix Mj=(Mj0, Mj1).
      // For the left half (generators), this is comparatively easy. The wavelet case
      // will cause much more effort due to the biorthogonalization equation
      //   Mj1 = (I-Mj0*(Mj0T^T))*Mj1c
      //       = Mj1c - Mj0*(Mj0T^T)*Mj1c
      
      InfiniteVector<double,IIndex> gcoeffs0, gcoeffs1, gcoeffs2;
      InfiniteVector<double,Index> psic; // psi_lambda^check
      
      const int ecode(lambda.e()[0]+2*lambda.e()[1]);
      
      switch(ecode) {
      case 0:
 	// generator
 	break;
      case 1:
 	// (1,0)-wavelet
 	break;
      case 2:
 	// (0,1)-wavelet
	
 	// First treat the identity part of the biorthogonalization equation, i.e.,
 	// compute the generator coefficients of the initial stable completion
	//   Mj1c = (I-Mj0*<Theta_{j+1},Lambda_j^tilde>^T)*Mj1
	// So, the wavelets which use one of the nonvaninshing boundary generators
	// will be modified to be zero at the endpoints.
	
//  	switch(lambda.p()) {
// 	case 0:
// 	  // psic_lambda is a tensor product of a generator and a wavelet on patch 0
// 	  basis00().reconstruct_1(IIndex(lambda.j(), 0, lambda.k()[0], &basis00()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis10().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis10()),
// 				  lambda.j()+1, gcoeffs1);

// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs1.begin()),
// 		   it2end(gcoeffs1.end());
// 		 it2 != it2end; ++it2) {
// 	      if (it2.index().k() == basis10().DeltaLmin()) {
// 		// interface generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					3,
// 					   typename Index::translation_type(it1.index().k(), 0),
// 					   this),
// 				     *it1 * *it2); // no factor!
// 	      } else {
// 		// patch generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   0,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
// 	      }
// 	    }
// 	  break;
// 	case 1:
// 	  // psic_lambda is a tensor product of a generator and a wavelet on patch 1
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 0, lambda.k()[0], &basis01()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis01()),
// 				  lambda.j()+1, gcoeffs1);

// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs1.begin()),
// 		   it2end(gcoeffs1.end());
// 		 it2 != it2end; ++it2) {
// 	      if (it2.index().k() == basis01().DeltaRmax(lambda.j()+1)) {
// 		// interface generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   3,
// 					   typename Index::translation_type(it1.index().k(), 0),
// 					   this),
// 				     *it1 * *it2); // no factor!
// 	      } else {
// 		// patch generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   1,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
// 	      }
// 	    }
// 	  break;
// 	case 2:
// 	  // psic_lambda is a tensor product of a generator and a wavelet on patch 2
// 	  basis10().reconstruct_1(IIndex(lambda.j(), 0, lambda.k()[0], &basis10()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis00().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis00()),
// 				  lambda.j()+1, gcoeffs1);

// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs1.begin()),
// 		   it2end(gcoeffs1.end());
// 		 it2 != it2end; ++it2) {
// 	      // (always!) a patch generator
// 	      psic.add_coefficient(Index(lambda.j()+1,
// 					 typename Index::type_type(),
// 					 2,
// 					 typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					 this),
// 				   *it1 * *it2);
// 	    }
// 	  break;
// 	case 4:
// 	  // psic_lambda decomposes into two tensor products of a generator and a wavelet,
// 	  // on patches 1 and 2
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 0, basis01().DeltaRmax(lambda.j()), &basis01()),
// 				  lambda.j()+1, gcoeffs0);
// 	  basis10().reconstruct_1(IIndex(lambda.j(), 0, basis10().DeltaLmin(), &basis10()),
// 				  lambda.j()+1, gcoeffs1);
// 	  basis01().reconstruct_1(IIndex(lambda.j(), 1, lambda.k()[1], &basis01()),
// 				  lambda.j()+1, gcoeffs2);
	  
// 	  // add the tensor product generators needed to construct psic_lambda
// 	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs0.begin()),
// 		 it1end(gcoeffs0.end());
// 	       it1 != it1end; ++it1)
// 	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs2.begin()),
// 		   it2end(gcoeffs2.end());
// 		 it2 != it2end; ++it2) {
// 	      if (it1.index().k() == basis01().DeltaRmax(lambda.j()+1)) {
// 		// interface generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   4,
// 					   typename Index::translation_type(0, it2.index().k()),
// 					   this),
// 				     *it1 * *it2); // no factor!
// 	      } else {
// 		// patch generator
// 		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   1,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
// 	      }
// 	    }
//  	  for (typename InfiniteVector<double,IIndex>::const_iterator it1(gcoeffs1.begin()),
//  		 it1end(gcoeffs1.end());
//  	       it1 != it1end; ++it1)
//  	    for (typename InfiniteVector<double,IIndex>::const_iterator it2(gcoeffs2.begin()),
//  		   it2end(gcoeffs2.end());
//  		 it2 != it2end; ++it2) {
//  	      if (it1.index().k() > basis10().DeltaLmin()) {
//  		// patch generator (interface generators already processed above)
//  		psic.add_coefficient(Index(lambda.j()+1,
// 					   typename Index::type_type(),
// 					   2,
// 					   typename Index::translation_type(it1.index().k(), it2.index().k()),
// 					   this),
// 				     *it1 * *it2);
//  	      }
//  	    }
// 	  break;
// 	} // end switch(lambda.p())
	
 	break;
      case 3:
 	// (1,1)-wavelet
 	break;
      }
      
      cout << "psic=" << endl << psic << endl;
      
      // compute c from psic
      c.swap(psic); // TODO: implement biorthogonalization
    }
  }
  
  template <class IBASIS>
  void
  LDomainBasis<IBASIS>::reconstruct_t_1(const Index& lambda,
					const int j,
					InfiniteVector<double, Index>& c) const {
  }

}

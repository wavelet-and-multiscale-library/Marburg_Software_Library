// implementation for ldomain_support.h

namespace WaveletTL
{
  template <class IBASIS>
  inline
  void support(const LDomainBasis<IBASIS>& basis,
	       const typename LDomainBasis<IBASIS>::Index& lambda,
	       typename LDomainBasis<IBASIS>::Support& supp)
  {
#if 1
    basis.support(lambda, supp);
#else
    // the following code has moved to LDomainBasis

    typedef typename LDomainBasis<IBASIS>::Index Index;

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

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[0],
		supp.xmax[0]);
	
	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[0],
		supp.ymax[0]);

	supp.xmin[1] = supp.xmin[2] = -1;

 	break;
      case 1:
 	// psi_lambda completely lives on patch 1

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[1],
		supp.xmax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[1],
		supp.ymax[1]);

	supp.xmin[0] = supp.xmin[2] = -1;

	break;
      case 2:
 	// psi_lambda completely lives on patch 2

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[2],
		supp.xmax[2]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[2],
		supp.ymax[2]);

	supp.xmin[0] = supp.xmin[1] = -1;

 	break;
      case 3:
  	// psi_lambda lives on patches 0 and 1
	
	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[0],
		supp.xmax[0]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       basis.basis1d().DeltaLmin(),
				       &basis.basis1d()),
		supp.ymin[0],
		supp.ymax[0]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[0],
				       &basis.basis1d()),
		supp.xmin[1],
		supp.xmax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       basis.basis1d().DeltaRmax(lambdaj),
				       &basis.basis1d()),
		supp.ymin[1],
		supp.ymax[1]);

	supp.xmin[2] = -1;

	break;
      case 4:
 	// psi_lambda lives on patches 1 and 2

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       basis.basis1d().DeltaRmax(lambdaj),
				       &basis.basis1d()),
		supp.xmin[1],
		supp.xmax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
		supp.ymin[1],
		supp.ymax[1]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       basis.basis1d().DeltaLmin(),
				       &basis.basis1d()),
		supp.xmin[2],
		supp.xmax[2]);

	support(basis.basis1d(),
		typename IBASIS::Index(lambdaj,
				       0,
				       lambda.k()[1],
				       &basis.basis1d()),
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
      InfiniteVector<double, Index> gcoeffs;
      basis.reconstruct_1(lambda, lambdaj+1, gcoeffs);
      
      typedef typename LDomainBasis<IBASIS>::Support Support;
      Support tempsupp;

      // initialize the support with an "empty" set
      for (int p = 0; p <= 2; p++) {
	supp.xmin[p] = -1;
      }

      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
	     itend(gcoeffs.end()); it != itend; ++it)
	{
	  // compute supp(psi_mu)
	  support(basis, it.index(), tempsupp);
	  
	  // for each patch p, update the corresponding support estimate
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
		supp.xmin[p] = std::min(supp.xmin[p], tempsupp.xmin[p]);
		supp.xmax[p] = std::max(supp.xmax[p], tempsupp.xmax[p]);
		supp.ymin[p] = std::min(supp.ymin[p], tempsupp.ymin[p]);
		supp.ymax[p] = std::max(supp.ymax[p], tempsupp.ymax[p]);
	      }
	    }
	  }
	}
    }
#endif
  }

  template <class IBASIS>
  bool intersect_supports(const LDomainBasis<IBASIS>& basis,
			  const typename LDomainBasis<IBASIS>::Index& lambda,
			  const typename LDomainBasis<IBASIS>::Support& supp2,
			  typename LDomainBasis<IBASIS>::Support& supp)
  {
    typedef typename LDomainBasis<IBASIS>::Support Support;

    Support supp1;
    support(basis, lambda, supp1);
    
    bool r = false;
    
    supp.j = std::max(supp1.j, supp2.j);
    const int diff1 = supp.j-supp1.j;
    const int diff2 = supp.j-supp2.j;
    
    for (int p = 0; p <= 2; p++) {
      if (supp1.xmin[p] != -1 && supp2.xmin[p] != -1) {
	// intersection of two nontrivial sets on patch p

	const int xmin1 = supp1.xmin[p] << diff1;
	const int xmax1 = supp1.xmax[p] << diff1;
	const int xmin2 = supp2.xmin[p] << diff2;
	const int xmax2 = supp2.xmax[p] << diff2;
	
	if (xmin1 < xmax2 && xmax1 > xmin2) {
	  // nontrivial intersection in x direction

	  const int ymin1 = supp1.ymin[p] << diff1;
	  const int ymax1 = supp1.ymax[p] << diff1;
	  const int ymin2 = supp2.ymin[p] << diff2;
	  const int ymax2 = supp2.ymax[p] << diff2;
	  
	  if (ymin1 < ymax2 && ymax1 > ymin2) {
	    // nontrivial intersection in y direction

	    supp.xmin[p] = std::max(xmin1, xmin2);
	    supp.xmax[p] = std::min(xmax1, xmax2);
	    supp.ymin[p] = std::max(ymin1, ymin2);
	    supp.ymax[p] = std::min(ymax1, ymax2);
	    
	    r = true;
	  } else {
	    supp.xmin[p] = -1;
	  }
	} else {
	  supp.xmin[p] = -1;
	}
      } else {
	supp.xmin[p] = -1;
      }
    }
    
    return r;
  }

  template <class IBASIS>
  bool intersect_supports(const LDomainBasis<IBASIS>& basis,
			  const typename LDomainBasis<IBASIS>::Index& lambda1,
			  const typename LDomainBasis<IBASIS>::Index& lambda2,
			  typename LDomainBasis<IBASIS>::Support& supp)
  {
    typedef typename LDomainBasis<IBASIS>::Support Support;

    Support supp2;
    support(basis, lambda2, supp2);

    return intersect_supports(basis, lambda1, supp2, supp);
  }
  
  template <class IBASIS>
  void intersecting_wavelets(const LDomainBasis<IBASIS>& basis,
			     const typename LDomainBasis<IBASIS>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename LDomainBasis<IBASIS>::Index>& intersecting)
  {
    typedef typename LDomainBasis<IBASIS>::Index Index;
    typedef typename LDomainBasis<IBASIS>::Support Support;

    Support supp, supp_lambda;
    support(basis, lambda, supp_lambda);
    
    // a brute force solution
    if (generators) {
      Index last_gen(last_generator<IBASIS>(&basis, j));
      for (Index mu = first_generator<IBASIS>(&basis, j);; ++mu) {
	if (intersect_supports(basis, mu, supp_lambda, supp))
	  intersecting.push_back(mu);
	if (mu == last_gen) break;
      }
    } else {
      Index last_wav(last_wavelet<IBASIS>(&basis, j));
      for (Index mu = first_wavelet<IBASIS>(&basis, j);; ++mu) {
	if (intersect_supports(basis, mu, supp_lambda, supp))
	  intersecting.push_back(mu);
	if (mu == last_wav) break;
      }
    }
  }
  
  template <class IBASIS>
  bool intersect_singular_support(const LDomainBasis<IBASIS>& basis,
				  const typename LDomainBasis<IBASIS>::Index& lambda,
				  const typename LDomainBasis<IBASIS>::Index& mu)
  {
    // we cheat a bit here: we return true if already the supports intersect (overestimate)
    return intersect_support(basis, lambda, mu);
  }

}

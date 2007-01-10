// implementation for ldomain_support.h

namespace WaveletTL
{
  template <class IBASIS>
  inline
  void support(const LDomainBasis<IBASIS>& basis,
	       const typename LDomainBasis<IBASIS>::Index& lambda,
	       typename LDomainBasis<IBASIS>::Support& supp)
  {
    basis.support(lambda, supp);
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
    
#if 1
    supp.j = supp1.j;
    int diff1 = 0;
    int diff2 = supp.j-supp2.j;
    if (supp2.j>supp1.j) {
      supp.j = supp2.j;
      diff1 = supp.j-supp1.j;
      diff2 = 0;
    }
#else
    // old version, slightly slower
    supp.j = std::max(supp1.j, supp2.j);
    const int diff1 = supp.j-supp1.j;
    const int diff2 = supp.j-supp2.j;
#endif
    
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

    intersecting.clear();

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
    typedef typename LDomainBasis<IBASIS>::Support Support;
    Support supp;
    return intersect_supports(basis, lambda, mu, supp);
  }

}

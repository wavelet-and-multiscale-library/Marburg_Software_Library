// implementation for ldomain_jl_support.h

namespace WaveletTL
{
  inline
  void support(const LDomainJLBasis& basis,
 	       const Index& lambda,
 	       Support& supp)
  {
    if (lambda.e()[0]+lambda.e()[1] == 0) {
      // generator
      supp.j = lambda.j();
      supp.xmin = lambda.k()[0]-1;
      supp.xmax = supp.xmin+2;
      supp.ymin = lambda.k()[1]-1;
      supp.ymax = supp.ymin+2;
    } else {
      // wavelet
      supp.j = lambda.j()+1;
      supp.xmin = 2*(lambda.k()[0]-1);
      supp.xmax = supp.xmin+4;
      supp.ymin = 2*(lambda.k()[1]-1);
      supp.ymax = supp.ymin+4;
    }
  }
 
  bool intersect_supports(const LDomainJLBasis& basis,
			  const Index& lambda,
			  const Support& supp2,
			  Support& supp)
  {
    Support supp1;
    support(basis, lambda, supp1);

    // determine support intersection granularity,
    // adjust single support granularities if necessary
    supp.j = std::max(supp1.j, supp2.j);
    const int diff1 = supp.j-supp1.j;
    const int diff2 = supp.j-supp2.j;

    const int xmin1 = supp1.xmin << diff1;
    const int xmax1 = supp1.xmax << diff1;
    const int xmin2 = supp2.xmin << diff2;
    const int xmax2 = supp2.xmax << diff2;

    if (xmin1 < xmax2 && xmax1 > xmin2) {
      // nontrivial intersection in x direction
      const int ymin1 = supp1.ymin << diff1;
      const int ymax1 = supp1.ymax << diff1;
      const int ymin2 = supp2.ymin << diff2;
      const int ymax2 = supp2.ymax << diff2;

      if (ymin1 < ymax2 && ymax1 > ymin2) {
	// nontrivial intersection in y direction
	supp.xmin = std::max(xmin1, xmin2);
	supp.xmax = std::min(xmax1, xmax2);
	supp.ymin = std::max(ymin1, ymin2);
	supp.ymax = std::min(ymax1, ymax2);
	return true;
      }
    }

    return false;
  }

  inline
  bool intersect_supports(const LDomainJLBasis& basis,
			  const Index& lambda,
			  const Index& mu,
			  Support& supp)
  {
    Support supp_mu;
    support(basis, mu, supp_mu);
    return intersect_supports(basis, lambda, supp_mu, supp);
  }

  void intersecting_wavelets(const LDomainJLBasis& basis,
			     const Index& lambda,
			     const int j, const bool generators,
			     std::list<Index>& intersecting) 
  {
    Support supp, supp_lambda;
    support(basis, lambda, supp_lambda);
    
#if 1
    // a brute force solution
    if (generators) {
      Index last_gen(last_generator(j));
      for (Index mu = first_generator(j);; ++mu) {
	if (intersect_supports(basis, mu, supp_lambda, supp))
	  intersecting.push_back(mu);
	if (mu == last_gen) break;
      }
    } else {
      Index last_wav(last_wavelet(j));
      for (Index mu = first_wavelet(j);; ++mu) {
	if (intersect_supports(basis, mu, supp_lambda, supp))
	  intersecting.push_back(mu);
	if (mu == last_wav) break;
      }
    }
#endif    
  }
 
  bool intersect_singular_support(const LDomainJLBasis& basis,
				  const Index& lambda,
				  const Index& mu)
  {
    // we cheat a bit here: we return true if already the supports intersect (overestimate)
    Support supp;
    return intersect_supports(basis, lambda, mu, supp);
  }
  

}

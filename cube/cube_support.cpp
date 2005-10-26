// implementation for cube_support.h

#include <utils/multiindex.h>

using MathTL::multi_degree;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM>
  void
  support(const CubeBasis<IBASIS,DIM>& basis,
	  const typename CubeBasis<IBASIS,DIM>::Index& lambda,
	  typename CubeBasis<IBASIS,DIM>::Support& supp)
  {
    const unsigned int jplus = multi_degree(lambda.e()) > 0 ? 1 : 0;
    supp.j = lambda.j() + jplus;
    for (unsigned int i(0); i < DIM; i++) {
      support(*basis.bases()[i],
	      typename IBASIS::Index(lambda.j(),
				     lambda.e()[i],
				     lambda.k()[i],
				     basis.bases()[i]),
	      supp.a[i], supp.b[i]);
      if (lambda.e()[i] == 0 && jplus > 0) {
	supp.a[i] *= 2;
	supp.b[i] *= 2;
      }
    }
  }

  template <class IBASIS, unsigned int DIM>
  bool
  intersect_supports(const CubeBasis<IBASIS,DIM>& basis,
		     const typename CubeBasis<IBASIS,DIM>::Index& lambda,
		     const typename CubeBasis<IBASIS,DIM>::Index& mu,
		     typename CubeBasis<IBASIS,DIM>::Support& supp)
  {
    typename CubeBasis<IBASIS,DIM>::Support supp_lambda;
    WaveletTL::support<IBASIS,DIM>(basis, lambda, supp_lambda);
      
    typename CubeBasis<IBASIS,DIM>::Support supp_mu;
    WaveletTL::support<IBASIS,DIM>(basis, mu, supp_mu);

    // determine support intersection granularity,
    // adjust single support granularities if necessary
    supp.j = std::max(supp_lambda.j, supp_mu.j);
    if (supp_lambda.j > supp_mu.j) {
      const int adjust = 1<<(supp_lambda.j-supp_mu.j);
      for (unsigned int i = 0; i < DIM; i++) {
	supp_mu.a[i] *= adjust;
	supp_mu.b[i] *= adjust;
      }
    } else {
      const int adjust = 1<<(supp_mu.j-supp_lambda.j);
      for (unsigned int i = 0; i < DIM; i++) {
	supp_lambda.a[i] *= adjust;
	supp_lambda.b[i] *= adjust;
      }
    }
    
    for (unsigned int i = 0; i < DIM; i++) {
      supp.a[i] = std::max(supp_lambda.a[i],supp_mu.a[i]);
      supp.b[i] = std::min(supp_lambda.b[i],supp_mu.b[i]);
      
      if (supp.a[i] >= supp.b[i])
	return false;
    }
    
    return true;
  }

  template <class IBASIS, unsigned int DIM>
  void intersecting_wavelets(const CubeBasis<IBASIS,DIM>& basis,
			     const typename CubeBasis<IBASIS,DIM>::Index& lambda,
			     const int j, const bool generators,
			     std::list<typename CubeBasis<IBASIS,DIM>::Index>& intersecting)
  {
    typedef typename CubeBasis<IBASIS,DIM>::Index Index;
    typedef typename CubeBasis<IBASIS,DIM>::Support Support;

    intersecting.clear();
      
    // a brute force solution
    Support supp;
    if (generators) {
      for (Index mu = first_generator<IBASIS,DIM>(&basis, j);; ++mu) {
	if (intersect_supports(basis, lambda, mu, supp))
	  intersecting.push_back(mu);
	if (mu == last_generator<IBASIS,DIM>(&basis, j)) break;
      }
    } else {
      for (Index mu = first_wavelet<IBASIS,DIM>(&basis, j);; ++mu) {
	if (intersect_supports(basis, lambda, mu, supp))
	  intersecting.push_back(mu);
	if (mu == last_wavelet<IBASIS,DIM>(&basis, j)) break;
      }
    }      
  }
    
  template <class IBASIS, unsigned int DIM>
  bool intersect_singular_support(const CubeBasis<IBASIS,DIM>& basis,
				  const typename CubeBasis<IBASIS,DIM>::Index& lambda,
				  const typename CubeBasis<IBASIS,DIM>::Index& nu)
  {
    return true; // TODO: implement this!
  }
}

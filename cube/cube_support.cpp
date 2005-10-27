// implementation for cube_support.h

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>

using MathTL::multi_degree;
using MathTL::FixedArray1D;

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

    intersecting.clear();

#if 1
    // the set of intersecting wavelets is a cartesian product from d sets from the 1D case,
    // so we only have to compute the relevant 1D indices
    typedef typename IBASIS::Index Index1D;
    FixedArray1D<std::list<Index1D>,DIM>
      intersecting_1d_generators, intersecting_1d_wavelets;

    for (unsigned int i = 0; i < DIM; i++) {
      intersecting_wavelets(*basis.bases()[i],
			    Index1D(lambda.j(),
				    lambda.e()[i],
				    lambda.k()[i],
				    basis.bases()[i]),
			    j, true, intersecting_1d_generators[i]);
      if (!(generators))
	intersecting_wavelets(*basis.bases()[i],
			      Index1D(lambda.j(),
				      lambda.e()[i],
				      lambda.k()[i],
				      basis.bases()[i]),
			      j, false, intersecting_1d_wavelets[i]);
    }

    // generate all relevant tensor product indices with either e=(0,...,0) or e!=(0,...,0)
    typedef std::list<FixedArray1D<Index1D,DIM> > list_type;
    list_type indices;
    FixedArray1D<Index1D,DIM> helpindex;
    if (DIM > 1 || (DIM == 1 && generators)) {
      for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[0].begin()),
	     itend(intersecting_1d_generators[0].end());
	   it != itend; ++it) {
	helpindex[0] = *it;
	indices.push_back(helpindex);
      }
    }
    if (!(generators)) {
      for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[0].begin()),
	     itend(intersecting_1d_wavelets[0].end());
	   it != itend; ++it) {
	helpindex[0] = *it;
	indices.push_back(helpindex);
      }
    }
    for (unsigned int i = 1; i < DIM; i++) {
      list_type sofar;
      sofar.swap(indices);
      for (typename list_type::const_iterator itready(sofar.begin()), itreadyend(sofar.end());
	   itready != itreadyend; ++itready) {
	helpindex = *itready;
	unsigned int esum = 0;
	for (unsigned int k = 0; k < i; k++)
	  esum += helpindex[i].e();
	if (generators || i < DIM-1 || (i == DIM-1 && esum > 0)) {
	  for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[i].begin()),
		 itend(intersecting_1d_generators[i].end());
	       it != itend; ++it) {
	    helpindex[i] = *it;
	    indices.push_back(helpindex);
	  }
	}
	if (!(generators)) {
	  for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[i].begin()),
		 itend(intersecting_1d_wavelets[i].end());
	       it != itend; ++it) {
	    helpindex[i] = *it;
	    indices.push_back(helpindex);
	  }
	}
      }
    }

    // compose the results
    typename Index::type_type help_e;
    typename Index::translation_type help_k;
    for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
	 it != itend; ++it) {
      for (unsigned int i = 0; i < DIM; i++) {
	help_e[i] = (*it)[i].e();
	help_k[i] = (*it)[i].k();
      }
      intersecting.push_back(Index(j, help_e, help_k, &basis));
    }
#else
    // a brute force solution
    typedef typename CubeBasis<IBASIS,DIM>::Support Support;
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
#endif      
  }
    
  template <class IBASIS, unsigned int DIM>
  bool intersect_singular_support(const CubeBasis<IBASIS,DIM>& basis,
				  const typename CubeBasis<IBASIS,DIM>::Index& lambda,
				  const typename CubeBasis<IBASIS,DIM>::Index& mu)
  {
    // we have intersection of the singular supports if and only if
    // one of the components have this property in one dimension
    typedef typename IBASIS::Index Index1D;
    for (unsigned int i = 0; i < DIM; i++) {
      if (intersect_singular_support
	  (*basis.bases()[i],
	   Index1D(lambda.j(), lambda.e()[i], lambda.k()[i], basis.bases()[i]),
	   Index1D(mu.j(), mu.e()[i], mu.k()[i], basis.bases()[i])))
	return true;
    }
    return false;
  }
}

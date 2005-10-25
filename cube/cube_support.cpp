// implementation for cube_support.h

#include <utils/multiindex.h>

using MathTL::multi_degree;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  void
  support(const CUBEBASIS& basis,
	  const typename CUBEBASIS::Index& lambda,
	  typename CUBEBASIS::Support& supp)
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

  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  bool
  intersect_supports_(const CUBEBASIS& basis,
		      const typename CUBEBASIS::Index& lambda,
		      const typename CUBEBASIS::Index& mu,
		      typename CUBEBASIS::Support& supp)
  {
    typename CUBEBASIS::Support supp_lambda;
    WaveletTL::support<IBASIS,DIM,CUBEBASIS>(basis, lambda, supp_lambda);

    typename CUBEBASIS::Support supp_mu;
    WaveletTL::support<IBASIS,DIM,CUBEBASIS>(basis, mu, supp_mu);


//     for (unsigned int i = 0; i< DIM; i++)
//       cout << supp_mu.j << " " << supp_mu.a[i] << " " << supp_mu.b[i] << endl;


    // determine integration granularity
    // and adjust granularities if necessary
    const int max_j = std::max(supp_lambda.j,supp_mu.j);
    if (supp_lambda.j > supp_mu.j) {
      for (unsigned int i = 0; i < DIM; i++) {      
	supp_mu.a[i] *= 1 << (supp_lambda.j-supp_mu.j);
	supp_mu.b[i] *= 1 << (supp_lambda.j-supp_mu.j);
      }
    }
    else if (supp_lambda.j < supp_mu.j){
      for (unsigned int i = 0; i < DIM; i++) {      
	supp_lambda.a[i] *= 1 << (supp_mu.j-supp_lambda.j);
	supp_lambda.b[i] *= 1 << (supp_mu.j-supp_lambda.j);
      }
    }
    
    supp.j = max_j;
    for (unsigned int i = 0; i < DIM; i++) {
      supp.a[i] = std::max(supp_lambda.a[i],supp_mu.a[i]);
      supp.b[i] = std::min(supp_lambda.b[i],supp_mu.b[i]);
      //cout << supp_intersect.a[i] << " " << supp_intersect.b[i] << endl;
      
      if (supp.a[i] >= supp.b[i])
	return false;
    }
    
    return true;
  
  }

}

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
}

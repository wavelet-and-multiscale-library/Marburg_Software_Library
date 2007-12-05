// implementation of cdf_basis.h

#include <cassert>
#include <numerics/cardinal_splines.h>

namespace WaveletTL
{
  template <int d, int dt>
  double
  CDFBasis<d,dt>::evaluate(const unsigned int derivative,
			   const RIndex& lambda,
			   const double x) const
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) // generator
      {
	r = (derivative == 0
	     ? MathTL::EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), x)
	     : MathTL::EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), x));
      }
    else // wavelet
      {
 	InfiniteVector<double, RIndex> gcoeffs;
	RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::reconstruct_1
	  (lambda, lambda.j()+1, gcoeffs);
 	for (typename InfiniteVector<double, RIndex>::const_iterator it(gcoeffs.begin());
 	     it != gcoeffs.end(); ++it)
 	  r += *it * evaluate(derivative, it.index(), x);
      }
    
    return r;
  }

}

// implementation for cdf_evaluate.h

#include <numerics/cardinal_splines.h>

namespace WaveletTL
{
  template <int d, int dT>
  double evaluate(const unsigned int derivative,
		  const RIndex& lambda,
		  const double x)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) // generator
      {
	r = (derivative == 0
	     ? EvaluateCardinalBSpline_td  (d, lambda.j(), lambda.k(), x)
	     : EvaluateCardinalBSpline_td_x(d, lambda.j(), lambda.k(), x));
      }
    else // wavelet
      {
 	InfiniteVector<double, RIndex> gcoeffs;
 	CDFBasis<d,dT> basis;
 	basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
 	for (typename InfiniteVector<double, RIndex>::const_iterator it(gcoeffs.begin());
 	     it != gcoeffs.end(); ++it)
 	  r += *it * evaluate<d,dT>(derivative, it.index(), x);
      }

    return r;
  }
}

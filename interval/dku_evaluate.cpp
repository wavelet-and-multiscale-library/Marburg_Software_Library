// implementation for dku_evaluate.h

#include <Rd/cdf_utils.h>
#include <algebra/infinite_vector.h>
#include <numerics/cardinal_splines.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <int d, int dT>
  double evaluate(const DKUBasis<d, dT>& basis, const unsigned int derivative,
		  const typename DKUBasis<d, dT>::Index& lambda,
		  const double x)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) // generator
      {
	if (lambda.k() <= basis.DeltaLmax())
	  {
 	    // left boundary generator
	    const int llklow  = 1-ell2<d>();
	    const int llkup   = basis.ellT_l()-1;

	    for (int i(0); i < llkup-llklow+1; i++)
	      {
 		double help(basis.CLA()(i, lambda.k()-basis.DeltaLmin()));
 		if (help != 0)
 		  r += help * (derivative == 0
			       ? EvaluateCardinalBSpline_td  (d, lambda.j(), llklow+i, x)
			       : EvaluateCardinalBSpline_td_x(d, lambda.j(), llklow+i, x));
	      }
	  }
	else
	  {
	    if (lambda.k() >= basis.DeltaRmin(lambda.j()))
	      {
		// right boundary generator
		const int rlklowh  = basis.ellT_r()-1;
		const int rlkuph   = 1-ell2<d>();
		
		for (int i(0); i < rlklowh-rlkuph+1; i++)
		  {
		    double help(basis.CRA()(i, basis.DeltaRmax(lambda.j())-lambda.k()));
		    if (help != 0)
		      r += help * (derivative == 0
				   ? EvaluateCardinalBSpline_td  (d, lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(rlkuph+i), x)
				   : EvaluateCardinalBSpline_td_x(d, lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(rlkuph+i), x));
		  }
	      }
	    else
	      {
		r = (derivative == 0
		     ? EvaluateCardinalBSpline_td  (d, lambda.j(), lambda.k(), x)
		     : EvaluateCardinalBSpline_td_x(d, lambda.j(), lambda.k(), x));
	      }
	  }
      }
    else // wavelet
      {
	InfiniteVector<double, typename DKUBasis<d, dT>::Index> gcoeffs;
	basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
	for (typename InfiniteVector<double, typename DKUBasis<d, dT>::Index>::const_iterator it(gcoeffs.begin());
	     it != gcoeffs.end(); ++it)
	  r += *it * evaluate(basis, derivative, it.index(), x);
      }

    return r;
  }
}

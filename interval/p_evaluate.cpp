// implementation for p_evaluate.h

#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <numerics/schoenberg_splines.h>

namespace WaveletTL
{
  template <int d, int dT>
  double evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
		  const typename PBasis<d,dT>::Index& lambda,
		  const double x)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) {
      // generator
      r = (derivative == 0
	   ? EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), x)
	   : EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), x));
    } else {
      // wavelet
      typedef typename PBasis<d,dT>::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
 	   it != gcoeffs.end(); ++it)
 	r += *it * evaluate(basis, derivative, it.index(), x);
    }
    
    return r;
  }
  
  template <int d, int dT>
  void
  evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
	   const typename PBasis<d,dT>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& values)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;

    if (lambda.e() == 0) {
      // generator
      for (unsigned int m(0); m < points.size(); m++)
	values[m] = (derivative == 0
		     ? EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), points[m])
		     : EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]));
    } else {
      // wavelet
      typedef typename PBasis<d,dT>::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      Array1D<double> help(points.size());
      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
	   it != gcoeffs.end(); ++it)
	{
	  evaluate(basis, derivative, it.index(), points, help);
	  for (unsigned int i = 0; i < points.size(); i++)
	    values[i] += *it * help[i];
	}
    }
  }

  template <int d, int dT>
  void evaluate(const PBasis<d,dT>& basis,
		const typename PBasis<d,dT>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    for (unsigned int i(0); i < npoints; i++) {
      funcvalues[i] = 0;
      dervalues[i] = 0;
    }

    if (lambda.e() == 0) {
      // generator
      for (unsigned int m(0); m < npoints; m++) {
	funcvalues[m] = EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), points[m]);
	dervalues[m]  = EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]);
      }
    } else {
      // wavelet
      typedef typename PBasis<d,dT>::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      Array1D<double> help1, help2;
      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
	   it != gcoeffs.end(); ++it)
	{
	  evaluate(basis, it.index(), points, help1, help2);
	  for (unsigned int i = 0; i < npoints; i++) {
	    funcvalues[i] += *it * help1[i];
	    dervalues[i]  += *it * help2[i];
	  }
	}
    }
  }
}

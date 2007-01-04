// implementation for ds_evaluate.h

#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <numerics/cardinal_splines.h>

namespace WaveletTL
{
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  double evaluate(const DSBasis<d,dT,BIO>& basis, const unsigned int derivative,
		  const typename DSBasis<d,dT,BIO>::Index& lambda,
		  const double x)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) {
      // generator
      const Matrix<double>& CLA(basis.get_CLA());
      if (lambda.k() < basis.DeltaLmin()+(int)CLA.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < CLA.row_dimension(); i++) {
	  double help(CLA(i, lambda.k()-basis.DeltaLmin()));
	  if (help != 0)
	    r += help * (derivative == 0
			 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, x)
			 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, x));
	}
      }	else {
	const Matrix<double>& CRA(basis.get_CRA());
	if (lambda.k() > basis.DeltaRmax(lambda.j())-(int)CRA.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < CRA.row_dimension(); i++) {
	    double help(CRA(i, basis.DeltaRmax(lambda.j())-lambda.k()));
	    if (help != 0)
	      r += help * (derivative == 0
			   ? EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x)
			   : EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x));
	  }
	} else {
	  // inner generator
	  r = (derivative == 0
	       ? EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), x)
	       : EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), x));
	}
      }
    } else {
      // wavelet
      typedef typename DSBasis<d,dT,BIO>::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
	   it != gcoeffs.end(); ++it)
	r += *it * evaluate(basis, derivative, it.index(), x);
    }
    
    return r;
  }
  
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  void
  evaluate(const DSBasis<d,dT,BIO>& basis, const unsigned int derivative,
	   const typename DSBasis<d,dT,BIO>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& values)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;

    if (lambda.e() == 0) {
      // generator
      const Matrix<double>& CLA(basis.get_CLA());
      if (lambda.k() < basis.DeltaLmin()+(int)CLA.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < CLA.row_dimension(); i++) {
	  const double help(CLA(i, lambda.k()-basis.DeltaLmin()));
	  // 	  if (help != 0)
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] += help * (derivative == 0
				 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m])
				 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, points[m]));
	}
      }	else {
	const Matrix<double>& CRA(basis.get_CRA());
	if (lambda.k() > basis.DeltaRmax(lambda.j())-(int)CRA.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < CRA.row_dimension(); i++) {
	    const double help(CRA(i, basis.DeltaRmax(lambda.j())-lambda.k()));
	    // 	    if (help != 0)
	    for (unsigned int m(0); m < points.size(); m++)
	      values[m] += help * (derivative == 0
				   ? EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m])
				   : EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m]));
	  }
	} else {
	  // inner generator
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = (derivative == 0
			 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), points[m])
			 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]));
	}
      }
    } else {
      // wavelet
      typedef typename DSBasis<d,dT,BIO>::Index Index;
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

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  void evaluate(const DSBasis<d,dT,BIO>& basis,
		const typename DSBasis<d,dT,BIO>::Index& lambda,
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
      const Matrix<double>& CLA(basis.get_CLA());
      if (lambda.k() < basis.DeltaLmin()+(int)CLA.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < CLA.row_dimension(); i++) {
	  const double help(CLA(i, lambda.k()-basis.DeltaLmin()));
	  // 	  if (help != 0)
	  for (unsigned int m(0); m < npoints; m++) {
	    funcvalues[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m]);
	    dervalues[m]  += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, points[m]);
	  }
	}
      }	else {
	const Matrix<double>& CRA(basis.get_CRA());
	if (lambda.k() > basis.DeltaRmax(lambda.j())-(int)CRA.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < CRA.row_dimension(); i++) {
	    const double help(CRA(i, basis.DeltaRmax(lambda.j())-lambda.k()));
	    // 	    if (help != 0)
	    for (unsigned int m(0); m < npoints; m++) {
	      funcvalues[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), points[m]);
	      dervalues[m]  += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), points[m]);
	    }
	  }
	} else {
	  // inner generator
	  for (unsigned int m(0); m < npoints; m++) {
	    funcvalues[m] = EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), points[m]);
	    dervalues[m]  = EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]);
	  }
	}
      }
    } else {
      // wavelet
      typedef typename DSBasis<d,dT,BIO>::Index Index;
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

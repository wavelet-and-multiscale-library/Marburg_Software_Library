// implementation for spline_evaluate.h

#include <interval/interval_bspline.h>

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor>
  SampledMapping<1>
  evaluate(const SplineBasis<d,dT,flavor>& basis,
	   const InfiniteVector<double, typename SplineBasis<d,dT,flavor>::Index>& coeffs,
	   const int resolution)
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization

    Grid<1> grid(0, 1, 1<<resolution);
    SampledMapping<1> result(grid); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it)
	jmax = std::max(it.index().j()+it.index().e(), jmax);

      // insert coefficients into a dense vector
      Vector<double> wcoeffs(basis.Deltasize(jmax));
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it) {
	// determine number of the wavelet
	typedef typename Vector<double>::size_type size_type;
	size_type number = 0;
	if (it.index().e() == 0) {
	  number = it.index().k()-basis.DeltaLmin();
	} else {
	  number = basis.Deltasize(it.index().j())+it.index().k()-basis.Nablamin();
	}
	wcoeffs[number] = *it;
      }
      
      // switch to generator representation
      Vector<double> gcoeffs(wcoeffs.size(), false);
      if (jmax == basis.j0())
	gcoeffs = wcoeffs;
      else
	basis.apply_Tj(jmax-1, wcoeffs, gcoeffs);
      
      Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++) {
	values[i] = 0;
	const double x = i*ldexp(1.0, -resolution);
	SchoenbergIntervalBSpline_td<d> sbs(jmax,0);
	for (unsigned int k = 0; k < gcoeffs.size(); k++) {
	  sbs.set_k(basis.DeltaLmin()+k);
	  values[i] += gcoeffs[k] * sbs.value(Point<1>(x));
	}
      }
      
      return SampledMapping<1>(grid, values);
    }
    
    return result;
  }



  template <int d, int dT, SplineBasisFlavor flavor>
  double evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
		  const typename SplineBasis<d,dT,flavor>::Index& lambda,
		  const double x)
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization
    assert(derivative <= 1); // we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	r = (derivative == 0
	     ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
							 (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							 1-x)
	     : -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
							  (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							  1-x));
      } else {
	r = (derivative == 0
	     ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), x)
	     : MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), x));
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      size_type number_lambda = basis.Deltasize(lambda.j())+lambda.k()-basis.Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      basis.apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	r += it->second * evaluate(basis, derivative, Index(lambda.j()+1, 0, basis.DeltaLmin()+it->first, &basis), x);
      }
    }
    
    return r;
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  void
  evaluate(const SplineBasis<d,dT,flavor>& basis, const unsigned int derivative,
	   const typename SplineBasis<d,dT,flavor>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& values)
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization
    assert(derivative <= 1); // we only support derivatives up to the first order

    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = (derivative == 0
		       ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								   (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								   1-points[m])
		       : -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
							    (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							    1-points[m]));
      } else {
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = (derivative == 0
		       ? MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								   lambda.k(),
								   points[m])
		       : MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								   lambda.k(),
								   points[m]));
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      size_type number_lambda = basis.Deltasize(lambda.j())+lambda.k()-basis.Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      basis.apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
      Array1D<double> help(points.size());
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(basis, derivative, Index(lambda.j()+1, 0, basis.DeltaLmin()+it->first, &basis), points, help);
	for (unsigned int i = 0; i < points.size(); i++)
	  values[i] += it->second * help[i];
      }
    }
  }
  
  template <int d, int dT, SplineBasisFlavor flavor>
  void evaluate(const SplineBasis<d,dT,flavor>& basis,
		const typename SplineBasis<d,dT,flavor>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    assert(flavor == P_construction); // the other cases are done in a template specialization

    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    for (unsigned int i(0); i < npoints; i++) {
      funcvalues[i] = 0;
      dervalues[i] = 0;
    }

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	for (unsigned int m(0); m < npoints; m++) {
	  funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								    (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								    1-points[m]);
	  dervalues[m]  = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								     (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								     1-points[m]);
	}
      } else {
	for (unsigned int m(0); m < npoints; m++) {
	  funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
								    lambda.k(),
								    points[m]);
	  dervalues[m]  = MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								    lambda.k(), 
								    points[m]);
	}
      }
    } else {
      // wavelet, switch to generator representation
      typedef typename Vector<double>::size_type size_type;
      size_type number_lambda = basis.Deltasize(lambda.j())+lambda.k()-basis.Nablamin();
      std::map<size_type,double> wc, gc;
      wc[number_lambda] = 1.0;
      basis.apply_Mj(lambda.j(), wc, gc);
      typedef typename SplineBasis<d,dT,flavor>::Index Index;
      Array1D<double> help1, help2;
      for (typename std::map<size_type,double>::const_iterator it(gc.begin());
	   it != gc.end(); ++it) {
	evaluate(basis, Index(lambda.j()+1, 0, basis.DeltaLmin()+it->first, &basis), points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += it->second * help1[i];
	  dervalues[i]  += it->second * help2[i];
	}
      }
    }
  }
  
}

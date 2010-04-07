// implementation for p_evaluate.h
#include <iostream>
#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <utils/array1d.h>
#include <numerics/schoenberg_splines.h>

namespace WaveletTL
{
  template <int d, int dT>
  SampledMapping<1>
  evaluate(const PBasis<d,dT>& basis,
	   const typename PBasis<d,dT>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    if (lambda.e() == 0) { // generator
      if (primal) {
	Grid<1> grid(0, 1, 1<<resolution);
 	MathTL::Array1D<double> values((1<<resolution)+1);
 	for (unsigned int i(0); i < values.size(); i++)
	  values[i] = evaluate(basis, 0, lambda, i*ldexp(1.0,-resolution));
	return SampledMapping<1>(grid, values);
      } else {
 	// dual
	int s0 = basis.get_s0();
	int s1 = basis.get_s1();
	const Matrix<double>& CLAT = basis.get_CLAT();
	// left boundary generator
	if (lambda.k() < basis.DeltaLTmin()+(int)CLAT.column_dimension()) {
	  if (s0 >= d-2) {
	    InfiniteVector<double, RIndex> coeffs;
	    for (unsigned int i(0); i < CLAT.row_dimension(); i++) {
	      double v(CLAT(i, lambda.k()-basis.DeltaLTmin()));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, 1-ell2T<d,dT>()+i), v);
	    }
	    return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	  }
	  else {
	    InfiniteVector<double, RIndex> coeffs;
	    for (unsigned int i(0); i < CLAT.row_dimension(); i++) {
	      double v(CLAT(i, lambda.k()-basis.DeltaLTmin()));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j()+1, 0, 1-ell2T<d,dT>()+i), v);
	    }
	    return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	  }
	}
	// no left boundary generator
	else {
	  const Matrix<double>& CRAT = basis.get_CRAT();
	  if (lambda.k() > basis.DeltaRTmax(lambda.j())-(int)CRAT.column_dimension()) {
	    if (s1 >= d-2) {
	      // right boundary generator
	      InfiniteVector<double, RIndex> coeffs;
	      for (unsigned int i(0); i < CRAT.row_dimension(); i++) {
		double v(CRAT(i, basis.DeltaRTmax(lambda.j())-lambda.k()));
		if (v != 0)
		  coeffs.set_coefficient(RIndex(lambda.j(), 0, (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2T<d,dT>()+i)), v);
	      }
	      return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	    }
	    else {
	      InfiniteVector<double, RIndex> coeffs;
	      for (unsigned int i(0); i < CRAT.row_dimension(); i++) {
		double v(CRAT(i, basis.DeltaRTmax(lambda.j())-lambda.k()));
		if (v != 0)
		  coeffs.set_coefficient(RIndex(lambda.j()+1, 0, (1<<(lambda.j()+1))-ell1<d>()-ell2<d>()-(1-ell2T<d,dT>()+i)), v);
	      }
	      return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	    }
	  }
	  // inner generator
	  else {
	    return basis.get_CDF_basis().evaluate(0, RIndex(lambda.j(), 0, lambda.k()),
						  primal, 0, 1, resolution);
	  }
	}
      }
    }
    else { // wavelet
      InfiniteVector<double, typename PBasis<d,dT>::Index> gcoeffs;
      if (primal)
	basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      else
	basis.reconstruct_t_1(lambda, lambda.j()+1, gcoeffs);
      
      return evaluate(basis, gcoeffs, primal, resolution);
    }
    
    return SampledMapping<1>(); // dummy return for the compiler
  }
  
  template <int d, int dT>
  SampledMapping<1>
  evaluate(const PBasis<d,dT>& basis,
	   const InfiniteVector<double, typename PBasis<d,dT>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    SampledMapping<1> result(Grid<1>(0, 1, 1<<resolution)); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename PBasis<d,dT>::Index Index;
      for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
	     itend(coeffs.end()); it != itend; ++it)
	jmax = std::max(it.index().j()+it.index().e(), jmax);

      // switch to generator representation
      InfiniteVector<double,Index> gcoeffs;
      if (primal)
	basis.reconstruct(coeffs,jmax,gcoeffs);
      else
	basis.reconstruct_t(coeffs,jmax,gcoeffs);

      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
	     itend(gcoeffs.end()); it != itend; ++it)
	result.add(*it, evaluate(basis, it.index(), primal, resolution));
    }
    
    return result;
  }
#include <iostream>
  template <int d, int dT>
  double evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
		  const typename PBasis<d,dT>::Index& lambda,
		  const double x)
  {
    assert(derivative <= 2); // we only support derivatives up to the second order
    double r = 0;
    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	switch (derivative){
	case 0: r= MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
							     (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							     1-x);
	  break;
	case 1: r=-MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
							     (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							     1-x); 
	  break;
	case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(),
							     (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							     1-x); 
	  break;
	}
      } else {
	switch (derivative){
	  case 0: r=MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), x);
	  break;
	  case 1: r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), x);
	  break;
	  case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(), lambda.k(), x);
	  break;
	}
      }
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

#ifdef _PP_AUSWERTUNG_DER_WAVELETS
  template <int d, int dT>
  Piecewise<double> expandAsPP(const PBasis<d,dT>& basis, const typename PBasis<d,dT>::Index& lambda)
  {
    assert(d <= 4); // we only support orders less then 4
    Piecewise<double> r;
    Polynomial<double> q;

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
	Polynomial<double> p;  // p(x) = 1-x
	p.set_coefficient(0, 1.0);
	p.set_coefficient(1, -1.0);
	r= MathTL::ExpandSchoenbergBspline<d>(lambda.j(),(1<<lambda.j())-d-lambda.k()-2*ell1<d>(),1);
	}
      else {
	r=MathTL::ExpandSchoenbergBspline<d>  (lambda.j(), lambda.k(),0);
	}
      }
    else {
      // wavelet
      typedef typename PBasis<d,dT>::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      for (typename InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
 	   it != gcoeffs.end(); ++it)
 	r += *it * expandAsPP(basis, it.index());
      }
    
    return r;

  }
#endif
  
  template <int d, int dT>
  void
  evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
	   const typename PBasis<d,dT>::Index& lambda,
	   const Array1D<double>& points, Array1D<double>& values)
  {   
    assert(derivative <= 2); // we only support derivatives up to the second order

#ifdef _PP_AUSWERTUNG_DER_WAVELETS // auswertung mit vorheriger umwandlung in PP
    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) 
	switch (derivative) {
	case 0:
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(),
								(1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								1-points[m]);
	  break;
	
	case 1: 
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								   (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								   1-points[m]);
	  break;

	case 2:
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(),
								   (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								   1-points[m]);
	  break;
	}
      else 
	switch (derivative) {
	  case 0: 
	    for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(),
								lambda.k(),
								points[m]);
	    break;
	
	case 1: 
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								  lambda.k(),
								  points[m]); 
	  break;

	case 2:
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(),
								  lambda.k(),
								  points[m]); 
	  break;

	}
      
    } else {
      // wavelet
      switch (derivative) {
	  case 0: 
	    for (unsigned int m(0); m < points.size(); m++){
	       values[m] = basis.wavelets[lambda.j()][lambda.k()](points[m]);
	     
	    }
	    break;
          case 1: 
	    for (unsigned int m(0); m < points.size(); m++){
	       values[m] = basis.wavelets[lambda.j()][lambda.k()].derivative(points[m]);
	     
	    }
	    break;
          case 2: 
	    for (unsigned int m(0); m < points.size(); m++){
	       values[m] = basis.wavelets[lambda.j()][lambda.k()].secondDerivative(points[m]);
	     
	    }
	    break;
          }
      }




#endif

#ifndef _PP_AUSWERTUNG_DER_WAVELETS  // auswertung ohne vorheriger umwandlung in PP
    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;

    if (lambda.e() == 0) {
      // generator
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) 
	switch (derivative) {
	case 0:
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(),
								(1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								1-points[m]);
	  break;
	
	case 1: 
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								   (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								   1-points[m]);
	  break;

	case 2:
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(),
								   (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
								   1-points[m]);
	  break;
	}
      else 
	switch (derivative) {
	  case 0: 
	    for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(),
								lambda.k(),
								points[m]);
	    break;
	
	case 1: 
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
								  lambda.k(),
								  points[m]); 
	  break;

	case 2:
	  for (unsigned int m(0); m < points.size(); m++)
	    values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(),
								  lambda.k(),
								  points[m]); 
	  break;

	}
      
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
#endif
  }

  template <int d, int dT>
  void evaluate(const PBasis<d,dT>& basis,
		const typename PBasis<d,dT>::Index& lambda,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
#ifdef _PP_AUSWERTUNG_DER_WAVELETS   // auswertung mit vorheriger umwandlung in PP
    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    //basis.setWavelets();
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
      // wavelet
	  for (unsigned int m(0); m < npoints; m++){
	     funcvalues[m] = basis.wavelets[lambda.j()][lambda.k()](points[m]);
	     dervalues[m] = basis.wavelets[lambda.j()][lambda.k()].derivative(points[m]);
	  }
      }
#endif   

#ifndef _PP_AUSWERTUNG_DER_WAVELETS   // auswertung ohne vorheriger umwandlung in PP
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
#endif
  }
}

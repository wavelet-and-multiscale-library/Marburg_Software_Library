// implementation for ds_evaluate.h

#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <numerics/cardinal_splines.h>

namespace WaveletTL
{
  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  SampledMapping<1>
  evaluate(const DSBasis<d,dT,BIO>& basis,
	   const typename DSBasis<d,dT,BIO>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    
    if (lambda.e() == 0) { // generator
      if (primal) {
	const Matrix<double>& CLA = basis.get_CLA();
	if (lambda.k() < basis.DeltaLmin()+(int)CLA.column_dimension()) {
	  // left boundary generator
	  InfiniteVector<double, RIndex> coeffs;
	  for(unsigned int i(0); i < CLA.row_dimension(); i++) {
	    double v(CLA(i, lambda.k()-basis.DeltaLmin()));
	    if (v != 0)
	      coeffs.set_coefficient(RIndex(lambda.j(), 0, 1-ell2<d>()+i), v);
	  }
	  return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	} else {
	  const Matrix<double>& CRA = basis.get_CRA();
	  if (lambda.k() > basis.DeltaRmax(lambda.j())-(int)CRA.column_dimension()) {
	    // right boundary generator
	    InfiniteVector<double, RIndex> coeffs;
	    for (unsigned int i(0); i < CRA.row_dimension(); i++) {
	      double v(CRA(i, basis.DeltaRmax(lambda.j())-lambda.k()));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i)), v);
	    }
	    return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	  } else {
	    // inner generator
	    return basis.get_CDF_basis().evaluate(0, RIndex(lambda.j(), 0, lambda.k()),
						  primal, 0, 1, resolution);
	  }
	}
      } else {
	// dual
	const Matrix<double>& CLAT = basis.get_CLAT();
	if (lambda.k() < basis.DeltaLTmin()+(int)CLAT.column_dimension()) {
	  // left boundary generator
	  InfiniteVector<double, RIndex> coeffs;
	  for (unsigned int i(0); i < CLAT.row_dimension(); i++) {
	    double v(CLAT(i, lambda.k()-basis.DeltaLTmin()));
	    if (v != 0)
	      coeffs.set_coefficient(RIndex(lambda.j(), 0, 1-ell2T<d,dT>()+i), v);
	  }
	  return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	} else {
	  const Matrix<double>& CRAT = basis.get_CRAT();
	  if (lambda.k() > basis.DeltaRTmax(lambda.j())-(int)CRAT.column_dimension()) {
	    // right boundary generator
	    InfiniteVector<double, RIndex> coeffs;
	    for (unsigned int i(0); i < CRAT.row_dimension(); i++) {
	      double v(CRAT(i, basis.DeltaRTmax(lambda.j())-lambda.k()));
	      if (v != 0)
		coeffs.set_coefficient(RIndex(lambda.j(), 0, (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2T<d,dT>()+i)), v);
	    }
	    return basis.get_CDF_basis().evaluate(0, coeffs, primal, 0, 1, resolution);
	  } else {
	    // inner generator
	    return basis.get_CDF_basis().evaluate(0, RIndex(lambda.j(), 0, lambda.k()),
						  primal, 0, 1, resolution);
	  }
	}
      }
    } else { // wavelet
      InfiniteVector<double, typename DSBasis<d,dT,BIO>::Index> gcoeffs;
      if (primal)
	basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      else
	basis.reconstruct_t_1(lambda, lambda.j()+1, gcoeffs);

      return evaluate(basis, gcoeffs, primal, resolution);
    }
    
    return SampledMapping<1>(); // dummy return for the compiler
  }

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  SampledMapping<1>
  evaluate(const DSBasis<d,dT,BIO>& basis,
	   const InfiniteVector<double, typename DSBasis<d,dT,BIO>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    
    SampledMapping<1> result(Grid<1>(0, 1, 1<<resolution)); // zero
   
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename DSBasis<d,dT,BIO>::Index Index;
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

  template <int d, int dT, DSBiorthogonalizationMethod BIO>
  double evaluate(const DSBasis<d,dT,BIO>& basis, const unsigned int derivative,
		  const typename DSBasis<d,dT,BIO>::Index& lambda,
		  const double x)
  {
   

    double r = 0;

    if (lambda.e() == 0) {
      // generator
      const Matrix<double>& CLA(basis.get_CLA());
      if (lambda.k() < basis.DeltaLmin()+(int)CLA.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < CLA.row_dimension(); i++) {
	  double help(CLA(i, lambda.k()-basis.DeltaLmin()));
	  if (help != 0) {
	    switch (derivative) {
	    case 0 :
	      r += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, x);
	      break;
	    case 1 :
	      r += help * EvaluateCardinalBSpline_td_x<d>  (lambda.j(), 1-ell2<d>()+i, x);
	      break;
	    case 2 :
	      r += help * EvaluateCardinalBSpline_td_x_sq<d>  (lambda.j(), 1-ell2<d>()+i, x);
	      break;
	    }// end switch
	  }
	    
// 	    r += help * (derivative == 0
// 			 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, x)
// 			 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, x));
	}
      }	else {
	const Matrix<double>& CRA(basis.get_CRA());
	if (lambda.k() > basis.DeltaRmax(lambda.j())-(int)CRA.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < CRA.row_dimension(); i++) {
	    double help(CRA(i, basis.DeltaRmax(lambda.j())-lambda.k()));
	    if (help != 0) {
	      switch (derivative) {
	      case 0 :
		r += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x);
		break;
	      case 1 :
		r += help * EvaluateCardinalBSpline_td_x<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x);
		break;
	      case 2 :
		r += help * EvaluateCardinalBSpline_td_x_sq<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x);
		break;
	      }// end switch
	      
	    }
	    // 	      r += help * (derivative == 0
// 			   ? EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x)
// 			   : EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-(d%2)-(1-ell2<d>()+i), x));
	  }
	} else {
	  // inner generator
	  switch (derivative) {
	  case 0 :
	    r = EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), x);
	    break;
	  case 1 :
	    r = EvaluateCardinalBSpline_td_x<d>  (lambda.j(), lambda.k(), x);
	    break;
	  case 2 :
	    r = EvaluateCardinalBSpline_td_x_sq<d>  (lambda.j(), lambda.k(), x);
	    break;
	  }// end switch

// 	  r = (derivative == 0
// 	       ? EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), x)
// 	       : EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), x));
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


    values.resize(points.size());
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;

    if (lambda.e() == 0) {
      // generator
      const Matrix<double>& CLA(basis.get_CLA());
      int z = 0;
      if (lambda.k() < basis.DeltaLmin()+(int)CLA.column_dimension()) {
	// left boundary generator
	for (unsigned int i(0); i < CLA.row_dimension(); i++) {
	  const double help(CLA(i, lambda.k()-basis.DeltaLmin()));
	  ++z;
	  //if (help != 0)
	  for (unsigned int m(0); m < points.size(); m++) {
	    //cout << "pints size = " << points.size() << endl;
	    switch (derivative) {
	    case 0 :
	      values[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m]);
	      break;
	    case 1 :
	      values[m] += help * EvaluateCardinalBSpline_td_x<d>  (lambda.j(), 1-ell2<d>()+i, points[m]);
	      break;
	    case 2 :
	      //	      cout << "z = " << z << endl;
// 	      cout << "p = " << points[m] << endl;
// 	      cout << "i " << i << " " << EvaluateCardinalBSpline_td_x_sq<d>  (lambda.j(), 1-ell2<d>()+i, points[m]) << endl;
  	      values[m] += help * EvaluateCardinalBSpline_td_x_sq<d>  (lambda.j(), 1-ell2<d>()+i, points[m]);
// 	      cout << help << " " << values[m] << endl;
	      //break;
	    };// end switch   
// 	    values[m] += help * (derivative == 0
// 				 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), 1-ell2<d>()+i, points[m])
// 				 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), 1-ell2<d>()+i, points[m]));
	  }
	}
      }	else {
	const Matrix<double>& CRA(basis.get_CRA());
	if (lambda.k() > basis.DeltaRmax(lambda.j())-(int)CRA.column_dimension()) {
	  // right boundary generator
	  for (unsigned int i(0); i < CRA.row_dimension(); i++) {
	    const double help(CRA(i, basis.DeltaRmax(lambda.j())-lambda.k()));
	    // 	    if (help != 0)
	    for (unsigned int m(0); m < points.size(); m++)
	      switch (derivative) {
	      case 0 :
		values[m] += help * EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m]);
		break;
	      case 1 :
		values[m] += help * EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m]);
		break;
	      case 2 :
		values[m] += help * EvaluateCardinalBSpline_td_x_sq<d>(lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m]);
		break;
	    }// end switch	    

// 	      values[m] += help * (derivative == 0
// 				   ? EvaluateCardinalBSpline_td<d>  (lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m])
// 				   : EvaluateCardinalBSpline_td_x<d>(lambda.j(), (1<<lambda.j())-ell1<d>()-ell2<d>()-(1-ell2<d>()+i), points[m]));
	  }
	} else {
	  // inner generator
	  for (unsigned int m(0); m < points.size(); m++)
	    switch (derivative) {
	    case 0 :
	      values[m] = EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), points[m]);
	      break;
	    case 1 :
	      values[m] = EvaluateCardinalBSpline_td_x<d>  (lambda.j(), lambda.k(), points[m]);
	      break;
	    case 2 :
	      values[m] = EvaluateCardinalBSpline_td_x_sq<d>  (lambda.j(), lambda.k(), points[m]);
	      break;
	    }// end switch	 
// 	    values[m] = (derivative == 0
// 			 ? EvaluateCardinalBSpline_td<d>  (lambda.j(), lambda.k(), points[m])
// 			 : EvaluateCardinalBSpline_td_x<d>(lambda.j(), lambda.k(), points[m]));
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

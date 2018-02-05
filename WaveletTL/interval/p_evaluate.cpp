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

#include "p_basis.h"
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
      InfiniteVector<double,int> gcoeffs;
      basis.reconstruct_1(lambda.j(), lambda.e(), lambda.k(), lambda.j()+1, gcoeffs); 
      const int Lmin(basis.DeltaLmin());
      for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
 	   it != gcoeffs.end(); ++it)
      {
          // gcoeffs contains only coeffs related to generators on level j+1
          // j_ = j+1; e_ = 0; k_ = DeltaLmin() + num
          r += *it * evaluate(basis, derivative, lambda.j()+1,0,Lmin+it.index(), x);
          //cout << "Index: " << it.index() << ", Wert: " << evaluate(basis, derivative, lambda.j()+1,0,Lmin+it.index(), x) << endl;
      }
    }
    return r;
  }
  
  template <int d, int dT>
  double evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
		  const int j, const int e, const int k,
		  const double x)
  {
    assert(derivative <= 2); // we only support derivatives up to the second order
    double r = 0;
    if (e == 0) {
      // generator
      if (k > (1<<j)-ell1<d>()-d) {
	switch (derivative){
	case 0: r= MathTL::EvaluateSchoenbergBSpline_td<d>  (j,
							     (1<<j)-d-k-2*ell1<d>(),
							     1-x);
	  break;
	case 1: r=-MathTL::EvaluateSchoenbergBSpline_td_x<d>(j,
							     (1<<j)-d-k-2*ell1<d>(),
							     1-x); 
	  break;
	case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(j,
							     (1<<j)-d-k-2*ell1<d>(),
							     1-x); 
	  break;
	}
      } else {
	switch (derivative){
	  case 0: r=MathTL::EvaluateSchoenbergBSpline_td<d>  (j, k, x);
	  break;
	  case 1: r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(j, k, x);
	  break;
	  case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(j, k, x);
	  break;
	}
      }
    } else {
      // wavelet
      InfiniteVector<double,int> gcoeffs;
      basis.reconstruct_1(j,e,k, j+1, gcoeffs); 
      const int Lmin(basis.DeltaLmin());
      for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
 	   it != gcoeffs.end(); ++it)
      {
          // gcoeffs contains only coeffs related to generators on level j+1
          // j_ = j+1; e_ = 0; k_ = DeltaLmin() + num
          r += *it * evaluate(basis, derivative, j+1,0,Lmin+it.index(), x);
      }
    }
 	//r += *it * evaluate(basis, derivative, it.index(), x);
    return r;
  }

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

    
    template <int d, int dT>
    void
    evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
            const typename PBasis<d,dT>::Index& lambda,
            const Array1D<double>& points, Array1D<double>& values)
    {   
        assert(derivative <= 2); // we only support derivatives up to the second order
        values.resize(points.size());
        for (unsigned int i(0); i < values.size(); i++)
            values[i] = 0;
        if (lambda.e() == 0) 
        {
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
        }
        else // wavelet
        {
            if(basis.get_evaluate_with_pre_computation()){ // with pre compuatation 
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
            else // not with pre computation
            {
                InfiniteVector<double,int> gcoeffs;
                basis.reconstruct_1(lambda.j(),lambda.e(), lambda.k(), lambda.j()+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help(points.size());
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, derivative, lambda.j()+1,0,Lmin+it.index(), points, help);
                    for (unsigned int i = 0; i < points.size(); i++)
                        values[i] += *it * help[i];
                }
                
            }
        }
    }
             
    template <int d, int dT>
    void evaluate(const PBasis<d,dT>& basis, const unsigned int derivative,
            const int j_, const int e_, const int k_,
            const Array1D<double>& points, Array1D<double>& values)
    {   
        assert(derivative <= 2); // we only support derivatives up to the second order
        values.resize(points.size());
        for (unsigned int i(0); i < values.size(); i++)
            values[i] = 0;
        if (e_ == 0) 
        {
            // generator
            if (k_ > (1<<j_)-ell1<d>()-d) 
                switch (derivative) {
                    case 0:
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(j_,
                                    (1<<j_)-d-k_-2*ell1<d>(),
                                    1-points[m]);
                        break;
                    case 1: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(j_,
                                    (1<<j_)-d-k_-2*ell1<d>(),
                                    1-points[m]);
                        break;
                    case 2:
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(j_,
                                    (1<<j_)-d-k_-2*ell1<d>(),
                                    1-points[m]);
                        break;
                }
            else 
                switch (derivative) {
                    case 0: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(j_,
                                    k_,
                                    points[m]);
                        break;
                    case 1: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td_x<d>(j_,
                                    k_,
                                    points[m]); 
                        break;
                    case 2:
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(j_,
                                    k_,
                                    points[m]); 
                        break;
                }
        }
        else // wavelet
        {
            if(basis.get_evaluate_with_pre_computation()){ // with pre compuatation 
                switch (derivative) {
                    case 0: 
                        for (unsigned int m(0); m < points.size(); m++){
                            values[m] = basis.wavelets[j_][k_](points[m]);
                        }
                        break;
                    case 1: 
                        for (unsigned int m(0); m < points.size(); m++){
                            values[m] = basis.wavelets[j_][k_].derivative(points[m]);
                        }
                        break;
                    case 2: 
                        for (unsigned int m(0); m < points.size(); m++){
                            values[m] = basis.wavelets[j_][k_].secondDerivative(points[m]);
                        }
                        break;
                }
            }
            else // not with pre computation
            {
                InfiniteVector<double,int> gcoeffs;
                basis.reconstruct_1(j_,e_,k_, j_+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help;
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, derivative, j_+1,0,Lmin+it.index(), points, help);
                    for (unsigned int i = 0; i < points.size(); i++)
                        values[i] += *it * help[i];
                }
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
        } 
        else {
            
            // wavelet
            if(basis.get_evaluate_with_pre_computation()) { // with pre compuatation 
                for (unsigned int m(0); m < npoints; m++){
                    funcvalues[m] = basis.wavelets[lambda.j()][lambda.k()](points[m]);
                    dervalues[m] = basis.wavelets[lambda.j()][lambda.k()].derivative(points[m]);
                }
            }
            else {  //without pre computation
                for (unsigned int i(0); i < npoints; i++) {
                    funcvalues[i] = 0;
                    dervalues[i] = 0;
                }
                // wavelet
                InfiniteVector<double,int> gcoeffs;
                basis.reconstruct_1(lambda.j(), lambda.e(), lambda.k(), lambda.j()+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help1, help2;
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, lambda.j()+1,0,Lmin+it.index(), points, help1, help2);
                    for (unsigned int i = 0; i < points.size(); i++) {
                        funcvalues[i] += *it * help1[i];
                        dervalues[i]  += *it * help2[i];
                    }
                }
            }
        }
    }
    
    template <int d, int dT>
    void evaluate(const PBasis<d,dT>& basis,
            const int j_, const int e_, const int k_,
            const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
    {
        const unsigned int npoints(points.size());
        funcvalues.resize(npoints);
        dervalues.resize(npoints);
        //basis.setWavelets();
        if (e_ == 0) {
            // generator
            if (k_ > (1<<j_)-ell1<d>()-d) {
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (j_,
                            (1<<j_)-d-k_-2*ell1<d>(),
                            1-points[m]);
                    dervalues[m]  = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(j_,
                            (1<<j_)-d-k_-2*ell1<d>(),
                            1-points[m]);
                }
            } else {
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (j_,
                            k_,
                            points[m]);
                    dervalues[m]  = MathTL::EvaluateSchoenbergBSpline_td_x<d>(j_,
                            k_, 
                            points[m]);
                }
            }
        } 
        else {
            // wavelet
            if(basis.get_evaluate_with_pre_computation()) { // with pre compuatation 
                for (unsigned int m(0); m < npoints; m++){
                    funcvalues[m] = basis.wavelets[j_][k_](points[m]);
                    dervalues[m] = basis.wavelets[j_][k_].derivative(points[m]);
                }
            }
            else {  //without pre computation
                for (unsigned int i(0); i < npoints; i++) {
                    funcvalues[i] = 0;
                    dervalues[i] = 0;
                }
                InfiniteVector<double,int> gcoeffs;
                basis.reconstruct_1(j_, e_, k_, j_+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help1, help2;
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, j_+1,0,Lmin+it.index(), points, help1, help2);
                    for (unsigned int i = 0; i < points.size(); i++) {
                        funcvalues[i] += *it * help1[i];
                        dervalues[i]  += *it * help2[i];
                    }
                }
            }
        }
    }
    
    
    //! Christoph: tuned versions. ATTENTION: ONLY LINEAR (d=2,dt=2) PRIMBS BASIS WITH HOMOG. BC's!!

    double evaluate_primbs_22_bc11(const int j, const int e, const int k, const double x)
    {

      double scale = (1u << (j/2));  // the scaling factor
      if ((j%2) == 1)
      {
        //scale *= sqrt(2);
        scale *= 1.41421356237309504880L;  // i.e., scale *= sqrt(2);
      }

      if (e==1)  // wavelet
      {
        double x_scaled = (1u << (j+1)) * x;
        int interval_choice = static_cast<int>(x_scaled);

        if (k < ((int)((1u << j) - 1)))   // inner or left boundary wavelet
        {
          if (k)                   // inner wavelet
          {
            if ( (interval_choice == (2*k - 2)) ||  (interval_choice == (2*k - 1)) )
            {
              return (scale*( x_scaled - (2*k - 2) )*0.25);
            }
            else if (interval_choice == (2*k))
            {
              return (scale*( 0.5 - 2*(x_scaled - 2*k) ));
            }
            else if (interval_choice == (2*k + 1))
            {
              return (scale*( -1.5 + 2*(x_scaled - 2*k - 1) ));
            }
            else if ( (interval_choice == (2*k + 2)) ||  (interval_choice == (2*k + 3)) )
            {
              return (scale*( 0.5 - 0.25*(x_scaled - 2*k - 2 ) ));
            }
            else
            {
              return 0.;
            }

          }
          else                     // left boundary wavelet
          {
            switch (interval_choice)
            {
              case 0:
                return ( scale * ( -1.25 * x_scaled ) );
                break;
              case 1:
                return ( scale * ( -1.25 + 2.75 * (x_scaled-1) ) );
                break;
              case 2:
              case 3:
                return ( scale * ( 1.5 - x_scaled + 2 ) );
                break;
              case 4:
              case 5:
                return ( scale * ( -0.5 + 0.25 * (x_scaled - 4) ) );
                break;
              default:
                return 0.;
            }
          }
        }
        else   // right boundary wavelet
        {
          int dummy_scale = 1u << (j+1);
          interval_choice -= (dummy_scale - 6);
          switch (interval_choice)
          {
            case (0):
            case (1):
              return ( scale * ( -0.25 *  (x_scaled - dummy_scale + 6) ) );
              break;
            case (2):
            case (3):
              return ( scale * ( -0.5 + x_scaled - dummy_scale + 4 ) );
              break;
            case (4):
              return ( scale * ( 1.5 - 2.75 * ( x_scaled - dummy_scale + 2 ) ) );
              break;
            case (5):
              return ( scale * ( -1.25 + 1.25 * (x_scaled - dummy_scale + 1) ) );
              break;
            default:
              return 0.;
          }
        }


      }    // Wavelet END
      else     // generator
      {
        double x_scaled = (1u << j) * x;

        if ( (x_scaled > (k-1)) && (x_scaled < (k+1)) )  // x is contained in support
        {
          if (x_scaled  < k)
          {
            return scale*(x_scaled - k + 1);
          }
          else
          {
            return scale*(k + 1 - x_scaled);
          }
        }
        else
        {
          return 0.;
        }

      }   // generator END
    }


    double evaluate_primbs_22_bc11_v3(const int j, const int e, const int k, const double x)
    {

      double scale = (1u << (j/2));  // the scaling factor
      if ((j%2) == 1)
      {
        //scale *= sqrt(2);
        scale *= 1.41421356237309504880L;  // i.e., scale *= sqrt(2);
      }

      if (e==1)  // wavelet
      {
        double x_scaled = (1u << (j+1)) * x;

        if (k < ((int)((1u << j) - 1)))   // inner or left boundary wavelet
        {
          if (k)                   // inner wavelet
          {
            int interval_choice = static_cast<int>(x_scaled) - 2*k + 2;
            switch (interval_choice)
            {
              case 0:
              case 1:
                return (scale*( x_scaled - (2*k - 2) )*0.25);
                break;
              case 2:
                return (scale*( 0.5 - 2*(x_scaled - 2*k) ));
                break;
              case 3:
                return (scale*( -1.5 + 2*(x_scaled - 2*k - 1) ));
                break;
              case 4:
              case 5:
                return (scale*( 0.5 - 0.25*(x_scaled - 2*k - 2 ) ));
                break;
              default:
                return 0.;
            }
          }
          else                     // left boundary wavelet
          {
            int interval_choice = static_cast<int>(x_scaled);
            switch (interval_choice)
            {
              case 0:
                return ( scale * ( -1.25 * x_scaled ) );
                break;
              case 1:
                return ( scale * ( -1.25 + 2.75 * (x_scaled-1) ) );
                break;
              case 2:
              case 3:
                return ( scale * ( 1.5 - x_scaled + 2 ) );
                break;
              case 4:
              case 5:
                return ( scale * ( -0.5 + 0.25 * (x_scaled - 4) ) );
                break;
              default:
                return 0.;
            }
          }
        }
        else   // right boundary wavelet
        {
          int dummy_scale = 1u << (j+1);
          int interval_choice = static_cast<int>(x_scaled) - dummy_scale + 6;

          switch (interval_choice)
          {
            case (0):
            case (1):
              return ( scale * ( -0.25 *  (x_scaled - dummy_scale + 6) ) );
              break;
            case (2):
            case (3):
              return ( scale * ( -0.5 + x_scaled - dummy_scale + 4 ) );
              break;
            case (4):
              return ( scale * ( 1.5 - 2.75 * ( x_scaled - dummy_scale + 2 ) ) );
              break;
            case (5):
              return ( scale * ( -1.25 + 1.25 * (x_scaled - dummy_scale + 1) ) );
              break;
            default:
              return 0.;
          }
        }


      }    // Wavelet END
      else     // generator
      {
        double x_scaled = (1u << j) * x;
        int interval_choice = static_cast<int>(x_scaled) - k + 1;

        switch (interval_choice)
        {
          case 0:
            return scale*(x_scaled - k + 1);
            break;
          case 1:
            return scale*(k + 1 - x_scaled);
            break;
          default:
            return 0.;
        }

      }   // generator END
    }



    double evaluate_primbs_22_bc11_deriv(const int j, const int e, const int k, const double x)
    {

      double scale = (1u << (j/2));  // the scaling factor
      if ((j%2) == 1)
      {
        //scale *= sqrt(2);
        scale *= 1.41421356237309504880L;  // i.e., scale *= sqrt(2);
      }

      if (e==1)  // wavelet
      {
        double x_scaled = (1u << (j+1)) * x;

        if (k < (int)(((1u << j) - 1)))   // inner or left boundary wavelet
        {
          if (k)                   // inner wavelet
          {
            int interval_choice = static_cast<int>(x_scaled) - 2*k + 2;
            switch (interval_choice)
            {
              case 0:
              case 1:
                return (scale*( 1u << (j-1) ));
                break;
              case 2:
                return (-scale*( 1u << (j+2) ));
                break;
              case 3:
                return (scale*( 1u << (j+2) ));
                break;
              case 4:
              case 5:
                return (-scale*( 1u << (j-1) ));
                break;
              default:
                return 0.;
            }
          }
          else                     // left boundary wavelet
          {
            int interval_choice = static_cast<int>(x_scaled);
            switch (interval_choice)
            {
              case 0:
                return (-scale*( 5 * (1u << (j-1)) ));
                break;
              case 1:
                return (scale*( 11 * (1u << (j-1)) ));
                break;
              case 2:
              case 3:
                return (-scale*( (1u << (j+1)) ));
                break;
              case 4:
              case 5:
                return (scale*( (1u << (j-1)) ));
                break;
              default:
                return 0.;
            }
          }
        }
        else   // right boundary wavelet
        {
          int dummy_scale = 1u << (j+1);
          int interval_choice = static_cast<int>(x_scaled) - dummy_scale + 6;

          switch (interval_choice)
          {
            case (0):
            case (1):
              return (-scale*( 1u << (j-1) ));
              break;
            case (2):
            case (3):
              return (scale*( 1u << (j+1) ));
              break;
            case (4):
              return (-scale*( 11 * (1u << (j-1)) ));
              break;
            case (5):
              return (scale*( 5 * (1u << (j-1)) ));
              break;
            default:
              return 0.;
          }
        }


      }    // Wavelet END
      else     // generator
      {
        double x_scaled = (1u << j) * x;
        int interval_choice = static_cast<int>(x_scaled) - k + 1;

        switch (interval_choice)
        {
          case 0:
            return scale*(1u << j);
            break;
          case 1:
            return -scale*(1u << j);
            break;
          default:
            return 0.;
        }

      }   // generator END
    }


    void evaluate_primbs_22_bc11_v3(const int j, const int e, const int k, const Array1D<double>& points, Array1D<double>& values)
    {

        values.resize(points.size());
        for (unsigned int i(0); i < values.size(); i++)
        {
          values[i] = 0.0;
        }

        for (unsigned int i(0); i < points.size(); i++)
        {
          values[i] = evaluate_primbs_22_bc11_v3(j, e, k, points[i]);
        }


    }
}

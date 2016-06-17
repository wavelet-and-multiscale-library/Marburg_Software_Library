// implementation for pq_evaluate.h
#include <iostream>
#include <Rd/r_index.h>
#include <Rd/cdf_utils.h>
#include <utils/array1d.h>
#include <numerics/schoenberg_splines.h>

namespace WaveletTL
{
  template <int d, int dT>
  SampledMapping<1>
  evaluate(const PQFrame<d,dT>& basis,
	   const typename PQFrame<d,dT>::Index& lambda,
	   const bool primal,
	   const int resolution,
           const int derivative)
  {
    if (lambda.e() == 0) { // generator
      if (primal) {
	Grid<1> grid(0, 1, 1<<resolution);
 	MathTL::Array1D<double> values((1<<resolution)+1);
        for (unsigned int i(0); i < values.size(); i++)
	  values[i] = evaluate(basis, derivative, lambda, i*ldexp(1.0,-resolution));
          
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
      InfiniteVector<double, typename PQFrame<d,dT>::Index> gcoeffs;
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
  evaluate(const PQFrame<d,dT>& basis,
	   const InfiniteVector<double, typename PQFrame<d,dT>::Index>& coeffs,
	   const bool primal,
	   const int resolution,
           const int derivative)
  {
    SampledMapping<1> result(Grid<1>(0, 1, 1<<resolution)); // zero
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef typename PQFrame<d,dT>::Index Index;
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
	result.add(*it, evaluate(basis, it.index(), primal, resolution, derivative));
    }
    
    return result;
  }
#include <iostream>

#include "pq_frame.h"
  template <int d, int dT>
  double evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
		  const typename PQFrame<d,dT>::Index& lambda,
		  const double x)
  {      
    assert(derivative <= 2); // we only support derivatives up to the second order
    double r = 0;
    
    const int pfkt = (d+1)/2;
    const int dhalfsmall = d/2;
    const double pfktrez = 1./pfkt;
    const double y = (1<<lambda.j())*x-lambda.k();
    const int lefthelper = lambda.k()+pfkt-1;
    const int righthelper = (1<<lambda.j())+dhalfsmall -lambda.k();
    const double leftside = 1./(pfkt+lambda.k());
    const double rightside = 1./righthelper;
    //cout << boundsupplength;
    
    
    if (lambda.e() == 0) {
      // generator
        //Für folgende Fälle noch nicht die Ableitungen angepasst @PHK
      if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {//right boundary quarks
	switch (derivative){
	case 0: r= MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
							     (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
							     1-x)*pow((1<<lambda.j())*(1-x)*rightside, lambda.p());
        //Scheint zu gehen; nochmal nachvollziehen und vor allem mit d=3 testen @PHK HIER WEITERMACHEN; mit bisherigem Ansatz werden vanishing moments zerstört @PHK
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
      } 
      
      
      else if (lambda.k() < -ell1<d>()){//left boundary quarks
	switch (derivative){
	  case 0: r=MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), x)*pow((1<<lambda.j())*x*leftside, lambda.p());
	  break;
	  case 1: r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), x);
	  break;
	  case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(), lambda.k(), x);
	  break;
	}
      }
      else{//inner quarks
        switch (derivative){
	  case 0: r=MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(), lambda.k(), x)*pow(y*pfktrez,lambda.p());
	  break;
	  case 1: (lambda.p()==0 ? r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), x)
                        : r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(), lambda.k(), x)*pow(y*pfktrez ,lambda.p())
                          + MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(), lambda.k(), x)
                          *lambda.p()*pow(y*pfktrez ,lambda.p()-1)*(1<<lambda.j())*pfktrez);
	  break;
	  case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(), lambda.k(), x);
	  break;    
        }
      }
       
    }
//        
     else {
      // wavelet
      InfiniteVector<double,int> gcoeffs;
      basis.reconstruct_1(lambda.p(), lambda.j(), lambda.e(), lambda.k(), lambda.j()+1, gcoeffs); 
      const int Lmin(basis.DeltaLmin());
      for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
 	   it != gcoeffs.end(); ++it)
      {
          // gcoeffs contains only coeffs related to generators on level j+1
          // j_ = j+1; e_ = 0; k_ = DeltaLmin() + num
          r += *it * evaluate(basis, derivative, lambda.p(), lambda.j()+1,0,Lmin+it.index(), x);
          //cout << "Stelle: " << it.index() << " Wert: " << evaluate(basis, derivative, lambda.p(), lambda.j()+1,0,Lmin+it.index(), x) << endl;
      }
    }
    return r;
  }
  
  template <int d, int dT>
  double evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
		  const int p, const int j, const int e, const int k,
		  const double x)
  {
    assert(derivative <= 2); // we only support derivatives up to the second order
    double r = 0;
    if (e == 0) {
      // generator
        
      const int pfkt = (d+1)/2;
      const int dhalfsmall = d/2;
      const double pfktrez = 1./pfkt;
      const double y = (1<<j)*x-k;
      const int lefthelper = k+pfkt-1;
      const int righthelper = (1<<j)+dhalfsmall -k;
      const double leftside = 1./(pfkt+k);
      const double rightside = 1./righthelper;
        
      if (k > (1<<j)-ell1<d>()-d) { //right boundary quarks
	switch (derivative){
            //Für folgende Fälle noch nicht die Ableitungen angepasst @PHK
	case 0: r= MathTL::EvaluateSchoenbergBSpline_td<d>(j,(1<<j)-d-k-2*ell1<d>(),1-x)
                                                          *pow((1<<j)*(1-x)*rightside, p);
        //cout << "Fall 1" << endl;
	  break;
	case 1: (p==0 ? r=-MathTL::EvaluateSchoenbergBSpline_td_x<d>(j,(1<<j)-d-k-2*ell1<d>(),1-x)
                      : r=-MathTL::EvaluateSchoenbergBSpline_td_x<d>(j,(1<<j)-d-k-2*ell1<d>(),1-x)*pow((1<<j)*(1-x)*rightside, p)
                          -MathTL::EvaluateSchoenbergBSpline_td<d>(j,(1<<j)-d-k-2*ell1<d>(),1-x)*p*pow((1<<j)*(1-x)*rightside, p-1)*(1<<j)*rightside); 
	  break;
	case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(j,(1<<j)-d-k-2*ell1<d>(),1-x); 
	  break;
	}
      }   
      else if (k < -ell1<d>()){ //left boundary quarks
	switch (derivative){
	  case 0: r=MathTL::EvaluateSchoenbergBSpline_td<d>(j, k, x)*pow((1<<j)*x*leftside, p);
          //cout << "Fall 2" << endl;
	  break;
	  case 1: (p==0 ? r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(j, k, x)
                        : r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(j, k, x)*pow((1<<j)*x*leftside, p)
                          +MathTL::EvaluateSchoenbergBSpline_td<d>(j, k, x)*p*pow((1<<j)*x*leftside, p-1)*(1<<j)*leftside);
	  break;
	  case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(j, k, x);
	  break;
          
	}
      }
        
      else {//inner quarks
	switch (derivative){
	  case 0: r=MathTL::EvaluateSchoenbergBSpline_td<d>(j, k, x)*pow(y*pfktrez ,p);
          //cout << "Fall 3" << endl;
	  break;
	  case 1: (p==0 ? r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(j, k, x) 
                        : r=MathTL::EvaluateSchoenbergBSpline_td_x<d>(j, k, x)*pow(y*pfktrez ,p)
                            + MathTL::EvaluateSchoenbergBSpline_td<d>(j, k, x)*p*pow(y*pfktrez ,p-1)*(1<<j)*pfktrez);
	  break;
	  //case 2: r=MathTL::EvaluateSchoenbergBSpline_td_xx<d>(j, k, x);
	  //break;
	}
      }
      
      
      
    } else {
      // wavelet
      InfiniteVector<double,int> gcoeffs;
      basis.reconstruct_1(p,j,e,k, j+1, gcoeffs); 
      const int Lmin(basis.DeltaLmin());
      for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
 	   it != gcoeffs.end(); ++it)
      {
          // gcoeffs contains only coeffs related to generators on level j+1
          // j_ = j+1; e_ = 0; k_ = DeltaLmin() + num
          r += *it * evaluate(basis, derivative, p,j+1,0,Lmin+it.index(), x);
      }
    }
 	//r += *it * evaluate(basis, derivative, it.index(), x);
    return r;
  }

  template <int d, int dT>
  Piecewise<double> expandAsPP(const PQFrame<d,dT>& basis, const typename PQFrame<d,dT>::Index& lambda)
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
      typedef typename PQFrame<d,dT>::Index Index;
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
    evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
            const typename PQFrame<d,dT>::Index& lambda,
            const Array1D<double>& points, Array1D<double>& values)
    {   
        assert(derivative <= 2); // we only support derivatives up to the second order
        values.resize(points.size());
        for (unsigned int i(0); i < values.size(); i++)
            values[i] = 0;
        if (lambda.e() == 0) 
        {
            // generator
            
            const int pfkt = (d+1)/2;
            const int dhalfsmall = d/2;
            const double pfktrez = 1./pfkt;
            const int lefthelper = lambda.k()+pfkt-1;
            const int righthelper = (1<<lambda.j())+dhalfsmall -lambda.k();
            const double leftside = 1./(pfkt+lambda.k());
            const double rightside = 1./righthelper;
            
            if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) 
                switch (derivative) {
                    case 0:
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(),
                                    (1<<lambda.j())-d-lambda.k()-2*ell1<d>(), 1-points[m])
                                    *pow((1<<lambda.j())*(1-points[m])*rightside, lambda.p());
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
            else if(lambda.k() < -ell1<d>())
                switch (derivative) {
                    case 0: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(),
                                    lambda.k(), points[m])
                                    *pow((1<<lambda.j())*points[m]*leftside, lambda.p());
                        break;
                    case 1: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
                                    lambda.k(), points[m]); 
                        break;
                    case 2:
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(),
                                    lambda.k(),
                                    points[m]);
                }
            else 
                switch (derivative) {
                    case 0: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(lambda.j(),
                                    lambda.k(), points[m])
                                    *pow(((1<<lambda.j())*points[m]-lambda.k())*pfktrez,lambda.p());
                        break;
                    case 1: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
                                    lambda.k(), points[m]); 
                        break;
                    case 2:
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td_xx<d>(lambda.j(),
                                    lambda.k(), points[m]); 
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
                basis.reconstruct_1(lambda.p(), lambda.j(),lambda.e(), lambda.k(), lambda.j()+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help(points.size());
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, derivative, lambda.p(), lambda.j()+1, 0, Lmin+it.index(), points, help);
                    for (unsigned int i = 0; i < points.size(); i++)
                        values[i] += *it * help[i];
                }
                
            }
        }
    }
             
    template <int d, int dT>
    void evaluate(const PQFrame<d,dT>& basis, const unsigned int derivative,
            const int p_, const int j_, const int e_, const int k_,
            const Array1D<double>& points, Array1D<double>& values)
    {   
        assert(derivative <= 2); // we only support derivatives up to the second order
        values.resize(points.size());
        for (unsigned int i(0); i < values.size(); i++)
            values[i] = 0;
        if (e_ == 0) 
        {
            // generator
            
            const int pfkt = (d+1)/2;
            const int dhalfsmall = d/2;
            const double pfktrez = 1./pfkt;
            const int lefthelper = k_+pfkt-1;
            const int righthelper = (1<<j_)+dhalfsmall -k_;
            const double leftside = 1./(pfkt+k_);
            const double rightside = 1./righthelper;
            
            if (k_ > (1<<j_)-ell1<d>()-d) 
                switch (derivative) {
                    case 0:
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(j_,
                                    (1<<j_)-d-k_-2*ell1<d>(),
                                    1-points[m])*pow((1<<j_)*(1-points[m])*rightside, p_);
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
            else if(k_ < -ell1<d>())
                switch (derivative) {
                    case 0: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(j_,
                                    k_, points[m])
                                    *pow((1<<j_)*points[m]*leftside, p_);
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
            
            
            else 
                switch (derivative) {
                    case 0: 
                        for (unsigned int m(0); m < points.size(); m++)
                            values[m] = MathTL::EvaluateSchoenbergBSpline_td<d>(j_,
                                    k_, points[m])
                                    *pow(((1<<j_)*points[m]-k_)*pfktrez,p_);
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
                basis.reconstruct_1(p_, j_,e_,k_, j_+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help;
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, derivative, p_, j_+1,0,Lmin+it.index(), points, help);
                    for (unsigned int i = 0; i < points.size(); i++)
                        values[i] += *it * help[i];
                }
            }
        }
    }

    template <int d, int dT>
    void evaluate(const PQFrame<d,dT>& basis,
            const typename PQFrame<d,dT>::Index& lambda,
            const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
    {
        const unsigned int npoints(points.size());
        funcvalues.resize(npoints);
        dervalues.resize(npoints);
        //basis.setWavelets();
        if (lambda.e() == 0) {
            // generator
            const int pfkt = (d+1)/2;
            const int dhalfsmall = d/2;
            const double pfktrez = 1./pfkt;
            const int lefthelper = lambda.k()+pfkt-1;
            const int righthelper = (1<<lambda.j())+dhalfsmall -lambda.k();
            const double leftside = 1./(pfkt+lambda.k());
            const double rightside = 1./righthelper;
            
            if (lambda.k() > (1<<lambda.j())-ell1<d>()-d) {
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
                            (1<<lambda.j())-d-lambda.k()-2*ell1<d>(), 1-points[m])
                            *pow((1<<lambda.j())*(1-points[m])*rightside, lambda.p());
                    dervalues[m]  = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
                            (1<<lambda.j())-d-lambda.k()-2*ell1<d>(),
                            1-points[m]);
                }
            } 
             else if(lambda.k() < -ell1<d>()){
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
                            lambda.k(), points[m])
                            *pow((1<<lambda.j())*points[m]*leftside, lambda.p());
                    dervalues[m]  = MathTL::EvaluateSchoenbergBSpline_td_x<d>(lambda.j(),
                            lambda.k(), 
                            points[m]);
                }
             }
            
             else {
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (lambda.j(),
                            lambda.k(), points[m])
                            *pow(((1<<lambda.j())*points[m]-lambda.k())*pfktrez,lambda.p());
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
                basis.reconstruct_1(lambda.p(), lambda.j(), lambda.e(), lambda.k(), lambda.j()+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help1, help2;
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, lambda.p(), lambda.j()+1,0,Lmin+it.index(), points, help1, help2);
                    for (unsigned int i = 0; i < points.size(); i++) {
                        funcvalues[i] += *it * help1[i];
                        dervalues[i]  += *it * help2[i];
                    }
                }
            }
        }
    }
    
    template <int d, int dT>
    void evaluate(const PQFrame<d,dT>& basis,
            const int p_, const int j_, const int e_, const int k_,
            const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
    {
        const unsigned int npoints(points.size());
        funcvalues.resize(npoints);
        dervalues.resize(npoints);
        //basis.setWavelets();
        if (e_ == 0) {
            // generator
            
            const int pfkt = (d+1)/2;
            const int dhalfsmall = d/2;
            const double pfktrez = 1./pfkt;
            const int lefthelper = k_+pfkt-1;
            const int righthelper = (1<<j_)+dhalfsmall -k_;
            const double leftside = 1./(pfkt+k_);
            const double rightside = 1./righthelper;
            
            if (k_ > (1<<j_)-ell1<d>()-d) {
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (j_,
                            (1<<j_)-d-k_-2*ell1<d>(), 1-points[m])
                            *pow((1<<j_)*(1-points[m])*rightside, p_);
                    dervalues[m]  = -MathTL::EvaluateSchoenbergBSpline_td_x<d>(j_,
                            (1<<j_)-d-k_-2*ell1<d>(),
                            1-points[m]);
                }
            }
            else if(k_ < -ell1<d>()){
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (j_,
                            k_, points[m])
                            *pow((1<<j_)*points[m]*leftside, p_);
                    dervalues[m]  = MathTL::EvaluateSchoenbergBSpline_td_x<d>(j_,
                            k_, 
                            points[m]);
                }
            }
            
            else {
                for (unsigned int m(0); m < npoints; m++) {
                    funcvalues[m] = MathTL::EvaluateSchoenbergBSpline_td<d>  (j_,
                            k_, points[m])
                            *pow(((1<<j_)*points[m]-k_)*pfktrez,p_);
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
                basis.reconstruct_1(p_, j_, e_, k_, j_+1, gcoeffs);
                const int Lmin(basis.DeltaLmin());
                Array1D<double> help1, help2;
                for (typename InfiniteVector<double,int>::const_iterator it(gcoeffs.begin());
                        it != gcoeffs.end(); ++it)
                {
                    evaluate(basis, p_, j_+1,0,Lmin+it.index(), points, help1, help2);
                    for (unsigned int i = 0; i < points.size(); i++) {
                        funcvalues[i] += *it * help1[i];
                        dervalues[i]  += *it * help2[i];
                    }
                }
            }
        }
    }
}

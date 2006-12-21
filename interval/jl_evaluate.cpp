// implementation for jl_evaluate.h

#include <geometry/grid.h>
#include <utils/array1d.h>
#include <numerics/schoenberg_splines.h>
#include <numerics/bezier.h>

using MathTL::Grid;

namespace WaveletTL
{
  SampledMapping<1>
  evaluate(const JLBasis& basis,
	   const JLBasis::Index& lambda,
	   const int resolution)
  {
    if (lambda.e() == 0) {
      // generator
      Grid<1> grid(0, 1, 1<<resolution);
      MathTL::Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++)
 	values[i] = evaluate(basis, 0, lambda, i*ldexp(1.0,-resolution));
      return SampledMapping<1>(grid, values);
    } else {
      // wavelet
#if 0
      // old code, uses JLBasis::reconstruct_1()
      InfiniteVector<double, JLBasis::Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);   
      return evaluate(basis, gcoeffs, resolution);
#else
      // new code after evaluate() works without reconstruct_1()
      Grid<1> grid(0, 1, 1<<resolution);
      MathTL::Array1D<double> values((1<<resolution)+1);
      for (unsigned int i(0); i < values.size(); i++)
 	values[i] = evaluate(basis, 0, lambda, i*ldexp(1.0,-resolution));
      return SampledMapping<1>(grid, values);
#endif
    }
    
    return SampledMapping<1>(); // dummy return for the compiler
  }
  
  SampledMapping<1>
  evaluate(const JLBasis& basis,
 	   const InfiniteVector<double, JLBasis::Index>& coeffs,
 	   const int resolution)
  {
    SampledMapping<1> result(Grid<1>(0, 1, 1<<resolution)); // zero
   
    if (coeffs.size() > 0) {
      // determine maximal level
      int jmax(0);
      typedef JLBasis::Index Index;
      for (InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
  	     itend(coeffs.end()); it != itend; ++it)
  	jmax = std::max(it.index().j()+it.index().e(), jmax);

      // switch to generator representation
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct(coeffs,jmax,gcoeffs);
      for (InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
  	     itend(gcoeffs.end()); it != itend; ++it)
  	result.add(*it, evaluate(basis, it.index(), resolution));
    }
    
    return result;
  }

  inline
  double evaluate(const JLBasis& basis, const unsigned int derivative,
  		  const JLBasis::Index& lambda,
  		  const double x)
  {
    return evaluate(derivative, lambda.j(), lambda.e(), lambda.c(), lambda.k(), x);
  }

  double evaluate(const unsigned int derivative,
		  const int j, const int e, const int c, const int k,
 		  const double x)
  {
    assert(derivative <= 1); // currently we only support derivatives up to the first order

    double r = 0;

    if (e == 0) {
      // generator
      const double factor = (c==0 ? sqrt(35./26.) : sqrt(105./2.)); // 1/||phi_c||_2
      r = factor * (derivative == 0
		    ? MathTL::EvaluateHermiteSpline_td  (c, j, k, x)
		    : MathTL::EvaluateHermiteSpline_td_x(c, j, k, x));
    } else {
      // wavelet
      const double phi0factor = sqrt(35./26.);
      const double phi1factor = sqrt(105./2.);
      if (c == 0) {
 	// type psi_0
 	// psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)
	const double factor = sqrt(35./352.);
	
 	int m = 2*k-1; // m-2k=-1
 	// phi_0(2x+1)
 	r += (-2.0)*M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x+1)
	r += (-21.0)*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
	
 	// m=2k <-> m-2k=0
 	m++;
 	// phi_0(2x)
	r += 4.0*M_SQRT1_2 * factor/phi0factor
	  * evaluate(derivative, j+1, 0, 0, m, x);
	
 	// m=2k+1 <-> m-2k=1
 	m++;
 	// phi_0(2x-1)
 	r += (-2.0)*M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x-1)
	r += 21.0*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
      } else { // c == 1
 	// type psi_1
 	// psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)
	const double factor = sqrt(35./48.);
	
 	int m = 2*k-1; // m-2k=-1
 	// phi_0(2x+1)
 	r += M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x+1)
	r += 9.0*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
	
 	// m=2k <-> m-2k=0
 	m++;
 	// phi_1(2x)
 	r += 12.0*M_SQRT1_2 * factor/phi1factor
 	  * evaluate(derivative, j+1, 0, 1, m, x);
	
 	// m=2k+1 <-> m-2k=1
 	m++;
 	// phi_0(2x-1)
 	r += (-1.0)*M_SQRT1_2 * factor/phi0factor
 	  * evaluate(derivative, j+1, 0, 0, m, x);
	// phi_1(2x-1)
	r += 9.0*M_SQRT1_2 * factor/phi1factor
	  * evaluate(derivative, j+1, 0, 1, m, x);
      }
    }

    return r;
  }
  
 
  void
  evaluate(const JLBasis& basis, const unsigned int derivative,
  	   const JLBasis::Index& lambda,
  	   const Array1D<double>& points, Array1D<double>& values)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    values.resize(points.size());

    if (lambda.e() == 0) {
      // generator
      const double factor = (lambda.c()==0 ? sqrt(35./26.) : sqrt(105./2.)); // 1/||phi_c||_2
      switch(derivative) {
      case 0:
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = factor
	    * MathTL::EvaluateHermiteSpline_td  (lambda.c(), lambda.j(), lambda.k(), points[m]);
	break;
      case 1:
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = factor
	    * MathTL::EvaluateHermiteSpline_td_x(lambda.c(), lambda.j(), lambda.k(), points[m]);
	break;
      default:
	break;
      }
    } else {
      // wavelet
      for (unsigned int i(0); i < values.size(); i++)
	values[i] = 0;
      
      typedef JLBasis::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      Array1D<double> help(points.size());
      for (InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
  	   it != gcoeffs.end(); ++it)
  	{
  	  evaluate(basis, derivative, it.index(), points, help);
  	  for (unsigned int i = 0; i < points.size(); i++)
  	    values[i] += *it * help[i];
  	}
    }
  }
  
  inline
  void evaluate(const int j, const int c, const int k,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    const double factor = (c==0 ? sqrt(35./26.) : sqrt(105./2.)); // 1/||phi_c||_2
    for (unsigned int m(0); m < npoints; m++) {
      funcvalues[m] = factor * MathTL::EvaluateHermiteSpline_td  (c, j, k, points[m]);
      dervalues[m]  = factor * MathTL::EvaluateHermiteSpline_td_x(c, j, k, points[m]);
    }
  }
  
  void evaluate(const JLBasis& basis,
  		const JLBasis::Index& lambda,
  		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    const unsigned int npoints(points.size());
    funcvalues.resize(npoints);
    dervalues.resize(npoints);
    if (lambda.e() == 0) {
      // generator
      const double factor = (lambda.c()==0 ? sqrt(35./26.) : sqrt(105./2.)); // 1/||phi_c||_2
      for (unsigned int m(0); m < npoints; m++) {
	funcvalues[m] = factor * MathTL::EvaluateHermiteSpline_td  (lambda.c(), lambda.j(), lambda.k(), points[m]);
	dervalues[m]  = factor * MathTL::EvaluateHermiteSpline_td_x(lambda.c(), lambda.j(), lambda.k(), points[m]);
      }
    } else {
      // wavelet
      for (unsigned int i(0); i < npoints; i++) {
	funcvalues[i] = 0;
	dervalues[i] = 0;
      }

      const double phi0factor = sqrt(35./26.);
      const double phi1factor = sqrt(105./2.);
      
      Array1D<double> help1, help2;
      if (lambda.c() == 0) {
	// type psi_0
	// psi_0(x) = -2*phi_0(2*x+1)+4*phi_0(2*x)-2*phi_0(2*x-1)-21*phi_1(2*x+1)+21*phi_1(2*x-1)
	const double factor = sqrt(35./352.);
	
	int m = 2*lambda.k()-1; // m-2k=-1
	// phi_0(2x+1)
	evaluate(lambda.j()+1, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-2.0)*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += (-2.0)*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x+1)
	evaluate(lambda.j()+1, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-21.0)*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += (-21.0)*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
	
	// m=2k <-> m-2k=0
	m++;
	// phi_0(2x)
	evaluate(lambda.j()+1, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 4.0*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += 4.0*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	
	// m=2k+1 <-> m-2k=1
	m++;
	// phi_0(2x-1)
	evaluate(lambda.j()+1, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-2.0)*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += (-2.0)*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x-1)
	evaluate(lambda.j()+1, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 21.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 21.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
      } else { // lambda.c() == 1
	// type psi_1
	// psi_1(x) = phi_0(2*x+1)-phi_0(2*x-1)+ 9*phi_1(2*x+1)+12*phi_1(2*x)+ 9*phi_1(2*x-1)
	const double factor = sqrt(35./48.);
	
	int m = 2*lambda.k()-1; // m-2k=-1
	// phi_0(2x+1)
	evaluate(lambda.j()+1, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x+1)
	evaluate(lambda.j()+1, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 9.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 9.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
	
	// m=2k <-> m-2k=0
	m++;
	// phi_1(2x)
	evaluate(lambda.j()+1, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 12.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 12.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
	
	// m=2k+1 <-> m-2k=1
	m++;
	// phi_0(2x-1)
	evaluate(lambda.j()+1, 0, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += (-1.0)*M_SQRT1_2 * factor/phi0factor * help1[i];
	  dervalues[i]  += (-1.0)*M_SQRT1_2 * factor/phi0factor * help2[i];
	}
	// phi_1(2x-1)
	evaluate(lambda.j()+1, 1, m, points, help1, help2);
	for (unsigned int i = 0; i < npoints; i++) {
	  funcvalues[i] += 9.0*M_SQRT1_2 * factor/phi1factor * help1[i];
	  dervalues[i]  += 9.0*M_SQRT1_2 * factor/phi1factor * help2[i];
	}
      }
    }
  }

}

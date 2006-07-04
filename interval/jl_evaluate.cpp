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
	   const bool primal,
	   const int resolution)
  {
    if (lambda.e() == 0) { // generator
      assert(primal);
      if (primal) {
 	Grid<1> grid(0, 1, 1<<resolution);
  	MathTL::Array1D<double> values((1<<resolution)+1);
   	for (unsigned int i(0); i < values.size(); i++)
  	  values[i] = evaluate(basis, 0, lambda, i*ldexp(1.0,-resolution));
 	return SampledMapping<1>(grid, values);
      } else {
  	// dual (unknown at the moment)
      }
    } else { // wavelet
      InfiniteVector<double, JLBasis::Index> gcoeffs;
      assert(primal);
      if (primal)
   	basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
//       else
//   	basis.reconstruct_t_1(lambda, lambda.j()+1, gcoeffs);

      return evaluate(basis, gcoeffs, primal, resolution);
    }
    
    return SampledMapping<1>(); // dummy return for the compiler
  }

  SampledMapping<1>
  evaluate(const JLBasis& basis,
	   const InfiniteVector<double, JLBasis::Index>& coeffs,
	   const bool primal,
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
      assert(primal);
      if (primal)
 	basis.reconstruct(coeffs,jmax,gcoeffs);
//       else
// 	basis.reconstruct_t(coeffs,jmax,gcoeffs);

      for (InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin()),
 	     itend(gcoeffs.end()); it != itend; ++it)
 	result.add(*it, evaluate(basis, it.index(), primal, resolution));
    }
    
    return result;
  }

  double evaluate(const JLBasis& basis, const unsigned int derivative,
 		  const JLBasis::Index& lambda,
 		  const double x)
  {
    assert(derivative <= 1); // currently we only support derivatives up to the first order

    double r = 0;

    if (lambda.e() == 0) {
      // generator

      assert(lambda.k() >= basis.DeltaLmin());
      assert(lambda.k() <= basis.DeltaRmax(lambda.j()));

      // phi_{j,s0},...,phi_{j,2^j-s1}             <-> 2^{j/2}phi_0(2^j*x-k), k=s0,...,s^j-s1
      // phi_{j,2^j-s1+1},...,phi_{j,2^{j+1}-s1+1} <-> 2^{j/2}phi_1(2^j*x-k), k=0,...,2^j
      const int first_half = (1<<lambda.j())-basis.get_s1();
      if (lambda.k() <= first_half) {
	// use phi_0
	const double phi_0_factor = sqrt(5./24.);
	r = phi_0_factor
	  * (derivative == 0
	     ? MathTL::EvaluateHermiteSpline_td  (0, lambda.j(), lambda.k(), x)
	     : MathTL::EvaluateHermiteSpline_td_x(0, lambda.j(), lambda.k(), x));
      } else {
	// use phi_1
	const int lambdak = lambda.k()-first_half-1;
	const double phi_1_factor = sqrt(15./8.);
	r = phi_1_factor
	  * (derivative == 0
	     ? MathTL::EvaluateHermiteSpline_td  (1, lambda.j(), lambdak, x)
	     : MathTL::EvaluateHermiteSpline_td_x(1, lambda.j(), lambdak, x));
      }
    } else {
      // wavelet

      assert(lambda.k() >= basis.Nablamin());
      assert(lambda.k() <= basis.Nablamax(lambda.j()));

      typedef JLBasis::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      for (InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
  	   it != gcoeffs.end(); ++it)
  	r += *it * evaluate(basis, derivative, it.index(), x);
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
    for (unsigned int i(0); i < values.size(); i++)
      values[i] = 0;

    if (lambda.e() == 0) {
      // generator

      assert(lambda.k() >= basis.DeltaLmin());
      assert(lambda.k() <= basis.DeltaRmax(lambda.j()));

      // phi_{j,s0},...,phi_{j,2^j-s1}             <-> 2^{j/2}phi_0(2^j*x-k), k=s0,...,s^j-s1
      // phi_{j,2^j-s1+1},...,phi_{j,2^{j+1}-s1+1} <-> 2^{j/2}phi_1(2^j*x-k), k=0,...,2^j
      
      const int first_half = (1<<lambda.j())-basis.get_s1();
      if (lambda.k() <= first_half) {
	// use phi_0
	const double phi_0_factor = sqrt(5./24.);
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = phi_0_factor
	    * (derivative == 0
	       ? MathTL::EvaluateHermiteSpline_td  (0, lambda.j(), lambda.k(), points[m])
	       : MathTL::EvaluateHermiteSpline_td_x(0, lambda.j(), lambda.k(), points[m]));
      } else {
	// use phi_1
	const double phi_1_factor = sqrt(15./8.);
	for (unsigned int m(0); m < points.size(); m++)
	  values[m] = phi_1_factor
	    * (derivative == 0
	       ? MathTL::EvaluateHermiteSpline_td  (1, lambda.j(), lambda.k()-first_half-1, points[m])
	       : MathTL::EvaluateHermiteSpline_td_x(1, lambda.j(), lambda.k()-first_half-1, points[m]));
      }
    } else {
      // wavelet

      assert(lambda.k() >= basis.Nablamin());
      assert(lambda.k() <= basis.Nablamax(lambda.j()));

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
  
  void evaluate(const JLBasis& basis,
 		const JLBasis::Index& lambda,
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

      assert(lambda.k() >= basis.DeltaLmin());
      assert(lambda.k() <= basis.DeltaRmax(lambda.j()));

      // phi_{j,s0},...,phi_{j,2^j-s1}             <-> 2^{j/2}phi_0(2^j*x-k), k=s0,...,s^j-s1
      // phi_{j,2^j-s1+1},...,phi_{j,2^{j+1}-s1+1} <-> 2^{j/2}phi_1(2^j*x-k), k=0,...,2^j
      
      const int first_half = (1<<lambda.j())-basis.get_s1();
      if (lambda.k() <= first_half) {
	// use phi_0
	const double phi_0_factor = sqrt(5./24.);
	for (unsigned int m(0); m < npoints; m++) {
	  funcvalues[m] = phi_0_factor
	    * MathTL::EvaluateHermiteSpline_td  (0, lambda.j(), lambda.k(), points[m]);
	  dervalues[m]  = phi_0_factor
	    * MathTL::EvaluateHermiteSpline_td_x(0, lambda.j(), lambda.k(), points[m]);
	}
      } else {
	// use phi_1
	const double phi_1_factor = sqrt(15./8.);
	for (unsigned int m(0); m < npoints; m++) {
	  funcvalues[m] = phi_1_factor
	    * MathTL::EvaluateHermiteSpline_td  (1, lambda.j(), lambda.k()-first_half-1, points[m]);
	  dervalues[m]  = phi_1_factor
	    * MathTL::EvaluateHermiteSpline_td_x(1, lambda.j(), lambda.k()-first_half-1, points[m]);
	}
      }
    } else {
      // wavelet

      assert(lambda.k() >= basis.Nablamin());
      assert(lambda.k() <= basis.Nablamax(lambda.j()));

      typedef JLBasis::Index Index;
      InfiniteVector<double,Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      Array1D<double> help1, help2;
      for (InfiniteVector<double,Index>::const_iterator it(gcoeffs.begin());
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

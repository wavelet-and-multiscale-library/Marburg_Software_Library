// implementation for jl_evaluate.h

#include <geometry/grid.h>
#include <utils/array1d.h>
#include <numerics/schoenberg_splines.h>
#include <numerics/bezier.h>
#include <interval/jl_utils.h>

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
 
  void
  evaluate(const JLBasis& basis, const unsigned int derivative,
  	   const JLBasis::Index& lambda,
  	   const Array1D<double>& points, Array1D<double>& values)
  {
    assert(derivative <= 1); // we only support derivatives up to the first order

    Array1D<double> dummy;
    if (derivative == 0)
      evaluate(lambda.j(), lambda.e(), lambda.c(), lambda.k(), points, values, dummy);
    else
      evaluate(lambda.j(), lambda.e(), lambda.c(), lambda.k(), points, dummy, values);
  }
  
  inline
  void evaluate(const JLBasis& basis,
  		const JLBasis::Index& lambda,
  		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues)
  {
    evaluate(lambda.j(), lambda.e(), lambda.c(), lambda.k(),
	     points, funcvalues, dervalues);
  }

}

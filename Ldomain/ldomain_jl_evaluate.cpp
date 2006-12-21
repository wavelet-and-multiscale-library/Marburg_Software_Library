// implementation for ldomain_jl_evaluate.h

#include <utils/fixed_array1d.h>
#include <interval/jl_utils.h>

using MathTL::FixedArray1D;
using MathTL::Point;

namespace WaveletTL
{
  Array1D<SampledMapping<2> >
  evaluate(const LDomainJLBasis& basis,
	   const Index& lambda,
	   const int resolution)
  {
    Array1D<SampledMapping<2> > r(3);

    // supp(psi_{j,e,c,k}) = 2^{-j}[k1-1,k1+1]x[k2-1,k2+1] cap Omega
    
    Index::type_type zero;
    if (lambda.e() == zero) {
      FixedArray1D<Array1D<double>,2> values;
      values[0].resize((1<<resolution)+1); // values in x-direction
      values[1].resize((1<<resolution)+1); // values in y-direction
      
      // patch Omega_0 = [-1,0]x[0,1]
      for (int i = 0; i <= 1<<resolution; i++)
	values[0][i] = values[1][i] = 0;
      if (lambda.k()[0] <= 0 && lambda.k()[1] >= 0) {
	for (int i = 0; i <= 1<<resolution; i++)
	  values[0][i] = evaluate(0, lambda.j(), 0, lambda.c()[0], lambda.k()[0],
				  -1.0+i*ldexp(1.0,-resolution));
	for (int i = 0; i <= 1<<resolution; i++)
	  values[1][i] = evaluate(0, lambda.j(), 0, lambda.c()[1], lambda.k()[1],
				  i*ldexp(1.0,-resolution));
      }
      r[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
      
      // patch Omega_1 = [-1,0]x[-1,0]
      for (int i = 0; i <= 1<<resolution; i++)
	values[0][i] = values[1][i] = 0;
      if (lambda.k()[0] <= 0 && lambda.k()[1] <= 0) {
	for (int i = 0; i <= 1<<resolution; i++)
	  values[0][i] = evaluate(0, lambda.j(), 0, lambda.c()[0], lambda.k()[0],
				  -1.0+i*ldexp(1.0,-resolution));
	for (int i = 0; i <= 1<<resolution; i++)
	  values[1][i] = evaluate(0, lambda.j(), 0, lambda.c()[1], lambda.k()[1],
				  -1.0+i*ldexp(1.0,-resolution));
      }
      r[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
      
      // patch Omega_2 = [0,1]x[-1,0]
      for (int i = 0; i <= 1<<resolution; i++)
	values[0][i] = values[1][i] = 0;
      if (lambda.k()[0] >= 0 && lambda.k()[1] <= 0) {
	for (int i = 0; i <= 1<<resolution; i++)
	  values[0][i] = evaluate(0, lambda.j(), 0, lambda.c()[0], lambda.k()[0],
				  i*ldexp(1.0,-resolution));
	for (int i = 0; i <= 1<<resolution; i++)
	  values[1][i] = evaluate(0, lambda.j(), 0, lambda.c()[1], lambda.k()[1],
				  -1.0+i*ldexp(1.0,-resolution));
      }
      r[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
    } else {
      InfiniteVector<double, Index> gcoeffs;
      basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
      return evaluate(basis, gcoeffs, resolution);
    }
    
    return r;
  }

  Array1D<SampledMapping<2> >
  evaluate(const LDomainJLBasis& basis,
 	   const InfiniteVector<double,Index>& coeffs,
 	   const int resolution)
  {
    // first prepare a plot of the zero function
    FixedArray1D<Array1D<double>,2> values;
    values[0].resize((1<<resolution)+1);
    values[1].resize((1<<resolution)+1);
    for (int i = 0; i <= 1<<resolution; i++) {
      values[0][i] = values[1][i] = 0;
    }
    Array1D<SampledMapping<2> > result(3);
    result[0] = SampledMapping<2>(Point<2>(-1, 0), Point<2>(0,1), values);
    result[1] = SampledMapping<2>(Point<2>(-1,-1), Point<2>(0,0), values);
    result[2] = SampledMapping<2>(Point<2>( 0,-1), Point<2>(1,0), values);
          
    // add all plots of the single functions
    for (InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
 	   itend(coeffs.end()); it != itend; ++it) {
      Array1D<SampledMapping<2> > temp(evaluate(basis, it.index(), resolution));
      result[0].add(*it, temp[0]);
      result[1].add(*it, temp[1]);
      result[2].add(*it, temp[2]);
    }
    
    return result;
  }

}

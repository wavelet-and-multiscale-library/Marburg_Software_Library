// implementation for gramian.h

namespace WaveletTL
{
  template <class RBASIS>
  double
  PeriodicIntervalGramian<RBASIS>::a(const typename WaveletBasis::Index& lambda,
				     const typename WaveletBasis::Index& mu,
				     const unsigned int p) const
  {
    double r = 0;
    
    if (WaveletBasis::intersect_supports(lambda, mu))
      {
	r = 42;
//  	// Set up Gauss points and weights for a composite quadrature formula:
//  	const unsigned int N_Gauss = d;
//  	const double h = ldexp(1.0, -supp.j);
//  	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values;
//  	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
//  	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
//  	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
//  	// - compute point values of the integrands
//   	evaluate(0, lambda, gauss_points, func1values);
//  	evaluate(0, mu, gauss_points, func2values);
	
//  	// - add all integral shares
//  	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
//  	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
//  	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
// 	    r += func1values[id] * func2values[id] * gauss_weight;
//  	  }
      }

    return r;
  }
  
}

// implementation for periodic_gramian.h

namespace WaveletTL
{
  template <class RBASIS>
  PeriodicIntervalGramian<RBASIS>::PeriodicIntervalGramian
  (const PeriodicBasis<RBASIS>& basis,
   const InfiniteVector<double, typename PeriodicBasis<RBASIS>::Index>& y)
    : basis_(basis), y_(y),
      normA(1.0),
      normAinv(ldexp(1.0, 2*(RBASIS::primal_polynomial_degree()-1))) // lower bound from [Bittner]
  {
  }
  
  template <class RBASIS>
  double
  PeriodicIntervalGramian<RBASIS>::a(const typename WaveletBasis::Index& lambda,
				     const typename WaveletBasis::Index& mu,
				     const unsigned int p) const
  {
    double r = 0;
    
    if (WaveletBasis::intersect_supports(lambda, mu))
      {
	// first we determine the support over which to integrate
	int j, k1, k2, length;
	if (lambda.j()+lambda.e() >= mu.j()+mu.e()) {
	  j = lambda.j()+lambda.e();
	  basis_.support(lambda, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	} else {
	  j = mu.j()+mu.e();
	  basis_.support(mu, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	}
	length = (k2 > k1 ? k2-k1 : k2+(1<<j)-k1); // number of subintervals
	
	// setup Gauss points and weights for a composite quadrature formula:
 	const unsigned int N_Gauss = RBASIS::primal_polynomial_degree()+1;
	const double h = ldexp(1.0, -j);

	Array1D<double> gauss_points (N_Gauss*(length)), func1values, func2values;
	int k = k1;
	for (int patch = 0; patch < length; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
	  for (unsigned int n = 0; n < N_Gauss; n++)
	    gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;

  	// - compute point values of the integrands
   	basis_.evaluate(0, lambda, gauss_points, func1values);
  	basis_.evaluate(0, mu, gauss_points, func2values);
	
  	// - add all integral shares
	for (int patch = k1, id = 0; patch < k1+length; patch++)
	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
	    r += func1values[id] * func2values[id] * GaussWeights[N_Gauss-1][n] * h;
	  }
      }
    
    return r;
  }
  
}

// implementation for ring_laplacian.h

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1>
  RingLaplacian<d,dt,s0,s1>::RingLaplacian(const RingBasis<d,dt,s0,s1>& basis,
					   const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y)
  {
  }
  
  template <int d, int dt, int s0, int s1>
  inline
  double
  RingLaplacian<d,dt,s0,s1>::a(const Index& lambda,
			       const Index& mu) const
  {
    // Some notes concerning the strategy are in order.
    // We have to compute the integral
    //   I = int_R nabla psi_lambda(x) nabla psi_mu(x) dx
    //     = int_0^1 int_0^1 nabla(psi^0_lambda(phi,s)) M(phi,s) (nabla(psi^0_mu(phi,s)))^T |det Dkappa(phi,s)| dphi ds
    // with the matrix
    //   M(phi,s) = D(kappa^{-1})(kappa(phi,s)) (D(kappa^{-1})(kappa(phi,s)))^T.
    // In the case of the ring-shaped domain, we have the chart
    //   kappa(phi,s) = r(s)*(cos(2*pi*phi),sin(2*pi*phi)), r(s)=r0+s*(r1-r0),
    // so that
    //   M(phi,s) = diag(1/(4*pi^2*r(s)^2), 1/(r1-r0)^2).
    // Using the fact that the wavelets psi_lambda and psi_mu are tensor products,
    //   psi_lambda(phi,s) = f1(phi)*f2(s)
    //   psi_mu    (phi,s) = g1(phi)*g2(s),
    // some calculation (cf. handritten notes+Maple) yields the following
    // representation as sum of products of 1D integrals:
    //   I = 1/(4*pi^2) * I1 * I3
    //       +1/(r1-r0)^2 * I2 * I4
    //       +1/4 * I2 * I3
    //       -1/(2*(r1-r0)) * I2 * (I5+I6)
    // with
    //   I1 = int_0^1 f1'(phi) g1'(phi) dphi
    //   I2 = int_0^1 f1(phi) g1(phi) dphi
    //   I3 = int_0^1 f2(s) g2(s) r(s)^{-2} ds
    //   I4 = int_0^1 f2'(s) g2'(s) ds
    //   I5 = int_0^1 f2(s) g2'(s) r(s)^{-1} ds
    //   I6 = int_0^1 f2'(s) g2(s) r(s)^{-1} ds
    // Note that {I1,I2} and {I3,I4,I5,I6} can be computed simultaneously,
    // using the same point evaluations of the integrand.
    // We do this in two helper subroutines.

    double I1, I2, I3, I4, I5, I6;

    angular_integrals
      (typename WaveletBasis::Index0(lambda.j(), lambda.e()[0], lambda.k()[0]),
       typename WaveletBasis::Index0(mu.j(), mu.e()[0], mu.k()[0]),
       I1, I2);
    
    radial_integrals
      (typename WaveletBasis::Index1(lambda.j(), lambda.e()[1], lambda.k()[1]),
       typename WaveletBasis::Index1(mu.j(), mu.e()[1], mu.k()[1]),
       I3, I4, I5, I6);
    
    return
      I1*I3/(4*M_PI*M_PI)
      + I2*I4/((basis_.r1()-basis_.r0())*(basis_.r1()-basis_.r0()))
      + I2*I3/4
      - I2*(I5+I6)/(2*(basis_.r1()-basis_.r0()));
  }
   
  template <int d, int dt, int s0, int s1>
  void
  RingLaplacian<d,dt,s0,s1>::angular_integrals
  (const Index0& lambda,
   const Index0& mu,
   double& I1, double& I2) const
  {
    // I1 = int_0^1 f1'(phi) g1'(phi) dphi
    // I2 = int_0^1 f1(phi) g1(phi) dphi

    I1 = I2 = 0;

    // lookup 1D integrals in the cache

    typename One_D_IntegralCache0::iterator col_lb(i12_cache.lower_bound(lambda));
    typename One_D_IntegralCache0::iterator col_it(col_lb);
    if (col_lb == i12_cache.end() ||
 	i12_cache.key_comp()(lambda,col_lb->first))
      {
 	// insert a new column
 	typedef typename One_D_IntegralCache0::value_type value_type;
 	col_it = i12_cache.insert(col_lb, value_type(lambda, Column1D_0()));
      }
    
    Column1D_0& col(col_it->second);
    
    typename Column1D_0::iterator lb(col.lower_bound(mu));
    typename Column1D_0::iterator it(lb);
    if (lb == col.end() ||
 	col.key_comp()(mu, lb->first))
      {
 	// cache miss

	if (basis_.basis0().intersect_supports(lambda, mu))
	  {
	    // first we determine the support over which to integrate
	    int j, k1, k2, length;
	    if (lambda.j()+lambda.e() >= mu.j()+mu.e()) {
	      j = lambda.j()+lambda.e();
	      basis_.basis0().support(lambda, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	    } else {
	      j = mu.j()+mu.e();
	      basis_.basis0().support(mu, k1, k2); // note: k2 may be less or equal to k1 in the case of overlap
	    }
	    length = (k2 > k1 ? k2-k1 : k2+(1<<j)-k1); // number of subintervals
	    
	    // setup Gauss points and weights for a composite quadrature formula:
	    const unsigned int N_Gauss = d+1;
	    const double h = 1.0/(1<<j);
	    
	    Array1D<double> gauss_points (N_Gauss*(length)),
	      func1values, func2values, der1values, der2values;
	    int k = k1;
	    for (int patch = 0; patch < length; patch++, k = dyadic_modulo(++k,j)) // work on 2^{-j}[k,k+1]
	      for (unsigned int n = 0; n < N_Gauss; n++)
		gauss_points[patch*N_Gauss+n] = h*(2*k+1+GaussPoints[N_Gauss-1][n])/2;
	    
	    // - compute point values of the integrands
	    basis_.basis0().evaluate(lambda, gauss_points, func1values, der1values);
	    basis_.basis0().evaluate(mu, gauss_points, func2values, der2values);
	    
	    // - add all integral shares
	    for (int patch = k1, id = 0; patch < k1+length; patch++)
	      for (unsigned int n = 0; n < N_Gauss; n++, id++) {
		const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
		I1 += gauss_weight * der1values[id] * der2values[id];
		I2 += gauss_weight * func1values[id] * func2values[id];
	      }
	  }

  	typedef typename Column1D_0::value_type value_type;
	FixedArray1D<double,2> help;
	help[0] = I1;
	help[1] = I2;
  	it = col.insert(lb, value_type(mu, help));
      }
    else
      {
 	// cache hit
 	I1 = it->second[0];
 	I2 = it->second[1];
      }
  }
  
  template <int d, int dt, int s0, int s1>
  void
  RingLaplacian<d,dt,s0,s1>::radial_integrals
  (const Index1& lambda,
   const Index1& mu,
   double& I3, double& I4, double& I5, double& I6) const
  {
    // I3 = int_0^1 f2(s) g2(s) r(s)^{-2} ds
    // I4 = int_0^1 f2'(s) g2'(s) ds
    // I5 = int_0^1 f2(s) g2'(s) r(s)^{-1} ds
    // I6 = int_0^1 f2'(s) g2(s) r(s)^{-1} ds

    I3 = I4 = I5 = I6 = 0;

    // lookup 1D integrals in the cache

    typename One_D_IntegralCache1::iterator col_lb(i3456_cache.lower_bound(lambda));
    typename One_D_IntegralCache1::iterator col_it(col_lb);
    if (col_lb == i3456_cache.end() ||
 	i3456_cache.key_comp()(lambda,col_lb->first))
      {
 	// insert a new column
 	typedef typename One_D_IntegralCache1::value_type value_type;
 	col_it = i3456_cache.insert(col_lb, value_type(lambda, Column1D_1()));
      }
    
    Column1D_1& col(col_it->second);
    
    typename Column1D_1::iterator lb(col.lower_bound(mu));
    typename Column1D_1::iterator it(lb);
    if (lb == col.end() ||
 	col.key_comp()(mu, lb->first))
      {
 	// cache miss

	// First we compute the support intersection of \psi_\lambda and \psi_\mu:
	typedef typename Basis1::Support Support;
	Support supp;
	
	if (basis_.basis1().intersect_supports(lambda, mu, supp))
	  {
	    // Set up Gauss points and weights for a composite quadrature formula:
	    const unsigned int N_Gauss = d+4;
	    const double h = 1.0/(1<<supp.j);
	    Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)),
	      func1values, func2values, der1values, der2values;
	    for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
	      for (unsigned int n = 0; n < N_Gauss; n++, id++)
		gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	    
	    // - compute point values of the integrands
	    basis_.basis1().evaluate(lambda, gauss_points, func1values, der1values);
	    basis_.basis1().evaluate(mu, gauss_points, func2values, der2values);
	    
	    // - add all integral shares
	    for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
	      for (unsigned int n = 0; n < N_Gauss; n++, id++) {
		const double rofs = basis_.r0()+gauss_points[id]*(basis_.r1()-basis_.r0());
		const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
		I3 += gauss_weight * func1values[id] * func2values[id] / (rofs*rofs);
		I4 += gauss_weight * der1values[id] * der2values[id];
		I5 += gauss_weight * func1values[id] * der2values[id] / rofs;
		I6 += gauss_weight * der1values[id] * func2values[id] / rofs;
	      }
	  }
	
  	typedef typename Column1D_1::value_type value_type;
	FixedArray1D<double,4> help;
	help[0] = I3;
	help[1] = I4;
	help[2] = I5;
	help[3] = I6;
  	it = col.insert(lb, value_type(mu, help));
      }
    else
      {
 	// cache hit
 	I3 = it->second[0];
 	I4 = it->second[1];
 	I5 = it->second[2];
 	I6 = it->second[3];
      }
  }
}

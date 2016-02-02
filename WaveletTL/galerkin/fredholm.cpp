// implementation for fredholm.h

namespace WaveletTL
{
  template <int d, int dT, int J0>
  FredholmIntegralOperator<d,dT,J0>::FredholmIntegralOperator
  (const WaveletBasis& basis,
   const InfiniteVector<double,Index>& y)
    : basis_(basis),
      y_(y),
      normA(0.0), normAinv(0.0)
  {
  }

  template <int d, int dT, int J0>
  FredholmIntegralOperator<d,dT,J0>::~FredholmIntegralOperator () {}

  template <int d, int dT, int J0>
  double
  FredholmIntegralOperator<d,dT,J0>::a
  (const Index& lambda,
   const Index& mu) const
  {
    double entry = 0;

    // compute the double integral
    //   <K^*Ku,v> = int_0^1 int_0^1 g(s,t) u(s) v(t) ds dt
    // with u=psi_lambda, v=psi_mu

    int j_lambda = lambda.j()+lambda.e();
    int j_mu = mu.j()+mu.e();
    int k1_lambda, k2_lambda, k1_mu, k2_mu;
    basis().support(lambda, k1_lambda, k2_lambda);
    basis().support(mu, k1_mu, k2_mu);

    // Set up Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = basis().primal_polynomial_degree()+1;
    const double h_lambda = ldexp(1.0, -j_lambda);
    const double h_mu = ldexp(1.0, -j_mu);
    Array1D<double> gp_lambda(N_Gauss*(k2_lambda-k1_lambda)),
      gp_mu(N_Gauss*(k2_mu-k1_mu)),
      psi_lambda_values, psi_mu_values;
    for (int patch = k1_lambda, id = 0; patch < k2_lambda; patch++)
      for (unsigned int n = 0; n < N_Gauss; n++, id++)
	gp_lambda[id] = h_lambda*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
    for (int patch = k1_mu, id = 0; patch < k2_mu; patch++)
      for (unsigned int n = 0; n < N_Gauss; n++, id++)
	gp_mu[id] = h_mu*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;

    // - compute point values of the integrands
    basis().evaluate(0, lambda, gp_lambda, psi_lambda_values);
    basis().evaluate(0, mu, gp_mu, psi_mu_values);

    // - add all integral shares
    for (int patch_lambda = k1_lambda, id2_lambda = 0; patch_lambda < k2_lambda; patch_lambda++)
      for (unsigned int n_lambda = 0; n_lambda < N_Gauss; n_lambda++, id2_lambda++) {
	const double gw_lambda = GaussWeights[N_Gauss-1][n_lambda] * h_lambda;
	for (int patch_mu = k1_mu, id2_mu = 0; patch_mu < k2_mu; patch_mu++)
	  for (unsigned int n_mu = 0; n_mu < N_Gauss; n_mu++, id2_mu++) {
	    const double gw_mu = GaussWeights[N_Gauss-1][n_mu] * h_mu;
	    entry +=
	      psi_lambda_values[id2_lambda] * gw_lambda
	      * psi_mu_values[id2_mu] * gw_mu
	      * g(gp_lambda[id2_lambda], gp_mu[id2_mu]);
	  }
      }
    
    return entry;
  }
  
  template <int d, int dT, int J0>
  VolterraIntegralOperator<d,dT,J0>::VolterraIntegralOperator
  (const WaveletBasis& basis,
   const InfiniteVector<double,Index>& y)
    : FredholmIntegralOperator<d,dT,J0>::FredholmIntegralOperator(basis, y)
  {
  }
}

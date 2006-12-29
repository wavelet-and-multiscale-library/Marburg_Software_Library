// implementation for ldomain_jl_laplacian.h

#include <Ldomain/ldomain_jl_support.h>
#include <numerics/gauss_data.h>
#include <interval/jl_utils.h>

namespace WaveletTL
{
  LDomainJLLaplacian::LDomainJLLaplacian(const WaveletBasis& basis,
					 const InfiniteVector<double,Index>& y)
    : basis_(basis), y_(y), normA(0.0), normAinv(0.0)
  {
  }
  
  double
  LDomainJLLaplacian::a(const Index& lambda,
			const Index& mu) const
  {
    double r = 0;
    
    // First we compute the support intersection of \psi_\lambda and \psi_\mu:
    typedef WaveletBasis::Support Support;
    Support supp;
    
    if (intersect_supports(basis_, lambda, mu, supp))
      {
	// use that both \psi_\lambda and \psi_\mu are a tensor product of 1D bases;
	// the entry of the Gramian on each square subpatch is a product of 2 1D integrals
	
	// iterate through the subsquares of supp and compute the integral shares
 	const double h = ldexp(1.0, -supp.j); // sidelength of the subsquare
	const int N_Gauss = 4;
	FixedArray1D<Array1D<double>,2> gauss_points, gauss_weights;
	for (int i = 0; i <= 1; i++) {
	  gauss_points[i].resize(N_Gauss);
	  gauss_weights[i].resize(N_Gauss);
	  for (int ii = 0; ii < N_Gauss; ii++)
	    gauss_weights[i][ii] = h*GaussWeights[N_Gauss-1][ii];
	}
//  	FixedArray1D<int,2> k;
//  	FixedArray1D<Array1D<double>,2>
//  	  psi_lambda_values,     // values of the components of psi_lambda at gauss_points[i]
//  	  psi_mu_values;         // -"-, for psi_mu
// 	Array1D<double> dummy;
//  	for (k[0] = supp.xmin; k[0] < supp.xmax; k[0]++) {
// 	  for (int ii = 0; ii < N_Gauss; ii++)
// 	    gauss_points[0][ii] = h*(2*k[0]+1+GaussPoints[N_Gauss-1][ii])/2.;
// 	  evaluate(lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0],
// 		   gauss_points[0], psi_lambda_values[0], dummy);
// 	  evaluate(mu.j(), mu.e()[0], mu.c()[0], mu.k()[0],
// 		   gauss_points[0], psi_mu_values[0], dummy);
// 	  double factor0 = 0;
// 	  for (int i0 = 0; i0 < N_Gauss; i0++)
// 	    factor0 += gauss_weights[0][i0] * psi_lambda_values[0][i0] * psi_mu_values[0][i0];
//  	  for (k[1] = supp.ymin; k[1] < supp.ymax; k[1]++) {
//  	    // check whether 2^{-supp.j}[k0,k0+1]x[k1,k1+1] is contained in Omega
//  	    if ((k[0] >= -(1<<supp.j) && k[0] < (1<<supp.j) && k[1] >= -(1<<supp.j) && k[1] < 0)
//  		|| (k[0] >= -(1<<supp.j) && k[0] < 0 && k[1] >= 0 && k[1] < (1<<supp.j))) {
// 	      for (int ii = 0; ii < N_Gauss; ii++)
// 		gauss_points[1][ii] = h*(2*k[1]+1+GaussPoints[N_Gauss-1][ii])/2.;
// 	      evaluate(lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1],
// 		       gauss_points[1], psi_lambda_values[1], dummy);
// 	      evaluate(mu.j(), mu.e()[1], mu.c()[1], mu.k()[1],
// 		       gauss_points[1], psi_mu_values[1], dummy);
// 	      double factor1 = 0;
// 	      for (int i1 = 0; i1 < N_Gauss; i1++)
// 		factor1 += gauss_weights[1][i1] * psi_lambda_values[1][i1] * psi_mu_values[1][i1];
	      
// 	      r += factor0 * factor1;
//  	    }
// 	  }
// 	}
      }
    
    return r;
  }
  
}

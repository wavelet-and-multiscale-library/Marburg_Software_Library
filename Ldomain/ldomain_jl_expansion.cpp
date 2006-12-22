// implementation for ldomain_jl_expansion.h

#include <Ldomain/ldomain_jl_support.h>
#include <numerics/gauss_data.h>
#include <interval/jl_utils.h>

namespace WaveletTL
{
  double
  integrate(const Function<2>* f,
	    const LDomainJLBasis& basis,
	    const Index& lambda)
  {
//     cout << "integrate() called with lambda=" << lambda << endl;

    double r = 0;
    
    typedef LDomainJLBasis::Support Support;
    Support supp;
    support(basis, lambda, supp);

    // iterate through the subsquares of supp and compute the integral shares
    const double h = ldexp(1.0, -supp.j); // sidelength of the subsquares
    const int N_Gauss = 6;
    FixedArray1D<Array1D<double>,2> gauss_points, gauss_weights;
    for (int i = 0; i <= 1; i++) {
      gauss_points[i].resize(N_Gauss);
      gauss_weights[i].resize(N_Gauss);
      for (int ii = 0; ii < N_Gauss; ii++)
	gauss_weights[i][ii] = h*GaussWeights[N_Gauss-1][ii];
    }
    FixedArray1D<int,2> k;
    FixedArray1D<Array1D<double>,2> psi_lambda_values; // values of the components of psi_lambda at gauss_points[i]
    for (k[0] = supp.xmin; k[0] < supp.xmax; k[0]++) {
      // evaluate first factor at Gauss points
      for (int ii = 0; ii < N_Gauss; ii++)
	gauss_points[0][ii] = h*(2*k[0]+1+GaussPoints[N_Gauss-1][ii])/2.;
      evaluate(0, lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0],
	       gauss_points[0], psi_lambda_values[0]);
      for (k[1] = supp.ymin; k[1] < supp.ymax; k[1]++) {
	// check whether 2^{-supp.j}[k0,k0+1]x[k1,k1+1] is contained in Omega
	if ((k[0] >= -(1<<supp.j) && k[0] < (1<<supp.j) && k[1] >= -(1<<supp.j) && k[1] < 0)
	    || (k[0] >= -(1<<supp.j) && k[0] < 0 && k[1] >= 0 && k[1] < (1<<supp.j))) {
//   	  cout << "in integrate(), [" << k[0]*h << "," << (k[0]+1)*h << "]x[" << k[1]*h << "," << (k[1]+1)*h << "] is in Omega" << endl;
	  // evaluate second factor at Gauss points
	  for (int ii = 0; ii < N_Gauss; ii++)
	    gauss_points[1][ii] = h*(2*k[1]+1+GaussPoints[N_Gauss-1][ii])/2.;
	  evaluate(0, lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1],
		   gauss_points[1], psi_lambda_values[1]);
//  	  cout << "gauss_points[0]=" << gauss_points[0] << endl;
// 	  cout << "psi_lambda_values[0]=" << psi_lambda_values[0] << endl;
//  	  cout << "gauss_points[1]=" << gauss_points[1] << endl;
// 	  cout << "psi_lambda_values[1]=" << psi_lambda_values[1] << endl;

	  // evaluate integral over the subsquare
	  Point<2> x;
	  for (int i0 = 0; i0 < N_Gauss; i0++) {
	    x[0] = gauss_points[0][i0];
	    for (int i1 = 0; i1 < N_Gauss; i1++) {
	      x[1] = gauss_points[1][i1];
	      r += f->value(x)
		* gauss_weights[0][i0] * psi_lambda_values[0][i0]
		* gauss_weights[1][i1] * psi_lambda_values[1][i1];
	    }
	  }
	} else {
// 	  cout << "in integrate(), [" << k[0]*h << "," << (k[0]+1)*h << "]x[" << k[1]*h << "," << (k[1]+1)*h << "] is NOT in Omega" << endl;
	}
      }
    }
      
    return r;
  }
  
  void
  expand(const Function<2>* f,
	 const LDomainJLBasis& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double,Index>& coeffs)
  {
    const int j0 = basis.j0();
    for (Index lambda = first_generator(j0);;++lambda)
      {
	const double coeff = integrate(f, basis, lambda);
	if (fabs(coeff) > 1e-14)
	  coeffs.set_coefficient(lambda, coeff);
 	if (lambda == last_wavelet(jmax))
 	  break;
      } 

    if (!primal) {
      // TODO!!!
    }
  }

}

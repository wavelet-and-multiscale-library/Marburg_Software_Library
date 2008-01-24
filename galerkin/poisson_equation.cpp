// implementation for poisson_equation.h

#include <numerics/quadrature.h>
#include <numerics/schoenberg_splines.h>
#include <interval/interval_bspline.h>
#include <numerics/eigenvalues.h>

using namespace MathTL;

namespace WaveletTL
{
  template <int d, int dT, int J0>
  PoissonEquation1D<d,dT,J0>::PoissonEquation1D
  (const WaveletBasis& basis,
   const InfiniteVector<double, typename WaveletBasis::Index>& y)
    : basis_(basis),
      A_(basis_, no_precond),
      y_(y),
      normA(0.0), normAinv(0.0)
  {
  }
  
  template <int d, int dT, int J0>
  inline
  double
  PoissonEquation1D<d,dT,J0>::a(const typename WaveletBasis::Index& lambda,
				const typename WaveletBasis::Index& nu) const
  {
    return a(lambda, nu, WaveletBasis::primal_polynomial_degree()*WaveletBasis::primal_polynomial_degree());
  }
  
  template <int d, int dT, int J0>
  double
  PoissonEquation1D<d,dT,J0>::a(const typename WaveletBasis::Index& lambda,
				const typename WaveletBasis::Index& nu,
				const unsigned int p) const
  {
    double r = 0;
    
    // first compute the support intersection of \psi_\lambda and \psi_\nu:
    typedef typename WaveletBasis::Support Support;
    Support supp;
    if (intersect_supports(basis_, lambda, nu, supp))
      {
#if 1
	// Set up Gauss points and weights for a composite quadrature formula:
	// (TODO: maybe use an instance of MathTL::QuadratureRule instead of computing
	// the Gauss points and weights)
	const unsigned int N_Gauss = (p+1)/2;
	const double h = ldexp(1.0, -supp.j);
	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), der1values, der2values;
	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
	// - compute point values of the integrands
	evaluate(basis_, 1, lambda, gauss_points, der1values);
	evaluate(basis_, 1, nu, gauss_points, der2values);

 	// - add all integral shares
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    r += der1values[id] * der2values[id] * gauss_weight;
 	  }
#else
 	// determine numbers of indices
	size_type number_lambda = 0, number_nu = 0;
	if (lambda.e() == 0) {
	  number_lambda = lambda.k()-basis_.DeltaLmin();
	} else {
	  number_lambda = basis_.Deltasize(lambda.j())+lambda.k()-basis_.Nablamin();
	}
	if (nu.e() == 0) {
	  number_nu = nu.k()-basis_.DeltaLmin();
	} else {
	  number_nu = basis_.Deltasize(nu.j())+nu.k()-basis_.Nablamin();
	}
	
	A_.set_level(std::max(lambda.j()+lambda.e(),nu.j()+nu.e()));
	return A_.get_entry(number_nu, number_lambda);
#endif
      }
    
    return r;
  }

  template <int d, int dT, int J0>
  double
  PoissonEquation1D<d,dT,J0>::norm_A() const
  {
    if (normA == 0.0) {
      FullLaplacian<d,dT,J0> A(basis(), energy);
      A.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }

    return normA;
  }
   
  template <int d, int dT, int J0>
  double
  PoissonEquation1D<d,dT,J0>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      FullLaplacian<d,dT,J0> A(basis(), energy);
      A.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(A, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }
    
    return normAinv;
  }
  
  template <int d, int dT, int J0>
  inline
  double
  PoissonEquation1D<d,dT,J0>::D(const typename WaveletBasis::Index& lambda) const
  {
#if 0
    // determine number of index lambda
    size_type number = 0;
    if (lambda.e() == 0) {
      number = lambda.k()-basis_.DeltaLmin();
    } else {
      number = basis_.Deltasize(lambda.j())+lambda.k()-basis_.Nablamin();
    }
    
    A_.set_level(lambda.j()+lambda.e());
    return sqrt(A_.diagonal(number));
#else
    return sqrt(a(lambda, lambda));
#endif
  }

}

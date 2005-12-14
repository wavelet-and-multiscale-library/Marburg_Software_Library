// implementation for gramian.h

#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/sparse_matrix.h>
#include <numerics/eigenvalues.h>
#include <numerics/gauss_data.h>
#include <numerics/iteratsolv.h>
#include <galerkin/galerkin_utils.h>

namespace WaveletTL
{
  template <class WBASIS>
  IntervalGramian<WBASIS>::IntervalGramian(const WBASIS& basis,
					   const InfiniteVector<double, typename WBASIS::Index>& y)
    : basis_(basis), y_(y), normA(0.0), normAinv(0.0)
  {
  }
  
  template <class WBASIS>
  double
  IntervalGramian<WBASIS>::a(const typename WaveletBasis::Index& lambda,
			     const typename WaveletBasis::Index& mu,
			     const unsigned int p) const
  {
    double r = 0;
    
    // First we compute the support intersection of \psi_\lambda and \psi_\mu:
    typedef typename WaveletBasis::Support Support;
    Support supp;
    
    if (intersect_supports(basis_, lambda, mu, supp))
      {
 	// Set up Gauss points and weights for a composite quadrature formula:
 	const unsigned int N_Gauss = basis_.dual_vanishing_moments();
 	const double h = ldexp(1.0, -supp.j);
 	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values;
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
 	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
 	// - compute point values of the integrands
  	evaluate(basis_, 0, lambda, gauss_points, func1values);
 	evaluate(basis_, 0, mu, gauss_points, func2values);
	
 	// - add all integral shares
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    r += func1values[id] * func2values[id] * gauss_weight;
 	  }
      }
    
    return r;
  }
  
  template <class WBASIS>
  double
  IntervalGramian<WBASIS>::norm_A() const
  {
    if (normA == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = 8;
      for (Index lambda = first_generator(&basis(), j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == last_wavelet(&basis(), jmax)) break;
      }
      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
#endif
    }
    
    return normA;
  }
  
  template <class WBASIS>
  double
  IntervalGramian<WBASIS>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      typedef typename WaveletBasis::Index Index;
      std::set<Index> Lambda;
      const int j0 = basis().j0();
      const int jmax = 8;
      for (Index lambda = first_generator(&basis(), j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == last_wavelet(&basis(), jmax)) break;
      }
      SparseMatrix<double> A_Lambda;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      
#if 1
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
#endif
    }
    
  }
}

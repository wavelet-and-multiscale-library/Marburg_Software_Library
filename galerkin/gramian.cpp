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
 	const unsigned int N_Gauss = basis_.primal_polynomial_degree()+1;
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

    return normAinv;
  }

  template <int d, int dT>
  IntervalGramian<SplineBasis<d,dT> >::IntervalGramian
  (const WaveletBasis& basis,
   const InfiniteVector<double, typename WaveletBasis::Index>& y)
    : basis_(basis),
      G_(basis_),
      y_(y),
      normA(0.0), normAinv(0.0)
  {
  }

  template <int d, int dT>
  double
  IntervalGramian<SplineBasis<d,dT> >::a
  (const typename WaveletBasis::Index& lambda,
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
	const unsigned int N_Gauss = p+1; //(p+1)/2;
	const double h = ldexp(1.0, -supp.j);
	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values;
	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
	// - compute point values of the integrands
	evaluate(basis_, 0, lambda, gauss_points, func1values);
	evaluate(basis_, 0, nu, gauss_points, func2values);

 	// - add all integral shares
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    r += func1values[id] * func2values[id] * gauss_weight;
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
	
	G_.set_level(std::max(lambda.j()+lambda.e(),nu.j()+nu.e()));
	return G_.get_entry(number_nu, number_lambda);
#endif
      }
    
    return r;
  }

  template <int d, int dT>
  double
  IntervalGramian<SplineBasis<d,dT> >::norm_A() const
  {
    if (normA == 0.0) {
      G_.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(G_, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }
    
    return normA;
  }
  
  template <int d, int dT>
  double
  IntervalGramian<SplineBasis<d,dT> >::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      G_.set_level(basis().j0()+4);
      double help;
      unsigned int iterations;
      LanczosIteration(G_, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
    }
    
    return normAinv;
  }
    
  template <int d, int dT>
  void
  IntervalGramian<SplineBasis<d,dT> >::add_level (const Index& lambda,
						  InfiniteVector<double, Index>& w, const int j,
						  const double factor,
						  const int J,
						  const CompressionStrategy strategy) const
  {
    // quick and dirty:
    // compute a full column of the stiffness matrix
    const int jmax = std::max(j+1, lambda.j()+lambda.e());
    G_.set_level(jmax);
    std::map<size_type,double> e_lambda, col_lambda;
    size_type number_lambda = lambda.k();
    if (lambda.e() == 0) {
      number_lambda -= basis_.DeltaLmin();
    } else {
      number_lambda += basis_.Deltasize(lambda.j())-basis_.Nablamin();
    }
    e_lambda[number_lambda] = 1.0;
    G_.apply(e_lambda, col_lambda);
    
    // extract the entries from level j
    if (j == basis_.j0()-1) {
      // "generator block"
      size_type startrow = 0;
      size_type endrow   = basis_.Deltasize(basis_.j0())-1;
      std::map<size_type,double>::const_iterator it(col_lambda.lower_bound(startrow));
      for (; it != col_lambda.end() && it->first <= endrow; ++it) {
	w.add_coefficient(Index(basis_.j0(), 0, basis_.DeltaLmin()+it->first, &basis_),
			  it->second * factor);
      }
    } else {
      // j>=j0, a "wavelet block"
      size_type startrow = basis_.Deltasize(j);
      size_type endrow   = basis_.Deltasize(j+1)-1;
      std::map<size_type,double>::const_iterator it(col_lambda.lower_bound(startrow));
      for (; it != col_lambda.end() && it->first <= endrow; ++it) {
	w.add_coefficient(Index(j, 1, basis_.Nablamin()+it->first-startrow, &basis_),
			  it->second * factor);
      }
    }
  }  
}

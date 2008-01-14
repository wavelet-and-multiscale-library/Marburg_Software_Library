// implementation for p_expansion.h

#include <set>
#include <list>

#include <utils/array1d.h>
#include <algebra/sparse_matrix.h>
#include <algebra/vector.h>
#include <numerics/gauss_data.h>
#include <numerics/iteratsolv.h>
#include <numerics/quadrature.h>
#include <galerkin/gramian.h>
#include <galerkin/cached_problem.h>
#include <interval/interval_bspline.h>
#include <adaptive/cdd1.h>

namespace WaveletTL
{
  template <int d, int dT>
  double integrate(const Function<1>* f,
		   const PBasis<d,dT>& basis,
		   const typename PBasis<d,dT>::Index& lambda)
  {
    double r = 0;
    
    const int j = lambda.j()+lambda.e();
    int k1, k2;
    support(basis, lambda, k1, k2);
    
    // setup Gauss points and weights for a composite quadrature formula:
    const unsigned int N_Gauss = 5;
    const double h = ldexp(1.0, -j);
    Array1D<double> gauss_points (N_Gauss*(k2-k1));
    for (int patch = k1; patch < k2; patch++) // refers to 2^{-j}[patch,patch+1]
      for (unsigned int n = 0; n < N_Gauss; n++)
	gauss_points[(patch-k1)*N_Gauss+n] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2;
    
    // add all integral shares
    for (unsigned int n = 0; n < N_Gauss; n++)
      {
	const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	for (int patch = k1; patch < k2; patch++)
	  {
	    const double t = gauss_points[(patch-k1)*N_Gauss+n];
	    
	    const double ft = f->value(Point<1>(t));
	    if (ft != 0)
	      r += ft
		* evaluate(basis, 0, lambda, t)
		* gauss_weight;
	  }
      }
    
    return r;
  }

  template <int d, int dT>
  double integrate(const PBasis<d,dT>& basis,
		   const typename PBasis<d,dT>::Index& lambda,
		   const typename PBasis<d,dT>::Index& mu)
  {
    double r = 0;
    
    // First we compute the support intersection of \psi_\lambda and \psi_\mu:
    typedef typename PBasis<d,dT>::Support Support;
    Support supp;

    if (intersect_supports(basis, lambda, mu, supp))
      {
 	// Set up Gauss points and weights for a composite quadrature formula:
 	const unsigned int N_Gauss = d;
 	const double h = ldexp(1.0, -supp.j);
 	Array1D<double> gauss_points (N_Gauss*(supp.k2-supp.k1)), func1values, func2values;
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++) // refers to 2^{-j}[patch,patch+1]
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++)
 	    gauss_points[id] = h*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
	
 	// - compute point values of the integrands
  	evaluate(basis, 0, lambda, gauss_points, func1values);
 	evaluate(basis, 0, mu, gauss_points, func2values);
	
 	// - add all integral shares
 	for (int patch = supp.k1, id = 0; patch < supp.k2; patch++)
 	  for (unsigned int n = 0; n < N_Gauss; n++, id++) {
 	    const double gauss_weight = GaussWeights[N_Gauss-1][n] * h;
	    r += func1values[id] * func2values[id] * gauss_weight;
 	  }
      }
    
    return r;
  }
  
  template <int d, int dT>
  void
  expand(const Function<1>* f,
	 const PBasis<d,dT>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename PBasis<d,dT>::Index>& coeffs)
  {
    typedef typename PBasis<d,dT>::Index Index;
    const int j0 = basis.j0();
    assert(jmax >= j0);
    
    coeffs.clear()

    for (Index lambda = first_generator(&basis, j0);;++lambda)
      {
 	coeffs.set_coefficient(lambda, integrate(f, basis, lambda));
 	if (lambda == last_wavelet(&basis, jmax))
 	  break;
      }

    if (!primal) {
#if 0
      IntervalGramian<PBasis<d,dT> > G(basis, coeffs);
      CachedProblem<IntervalGramian<PBasis<d,dT> > > GC(&G);
      InfiniteVector<double, typename PBasis<d,dT>::Index> x;
      CDD1_SOLVE(GC, 1e-6, x, jmax);
      coeffs.swap(x);
#else
      // setup active index set
      std::set<Index> Lambda;
      for (Index lambda = first_generator(&basis, j0);; ++lambda) {
 	Lambda.insert(lambda);
 	if (lambda == last_wavelet(&basis, jmax)) break;
      }
      
      // setup Gramian A_Lambda
      SparseMatrix<double> A_Lambda(Lambda.size());
      typedef typename SparseMatrix<double>::size_type size_type;     
      size_type row = 0;
      for (typename std::set<Index>::const_iterator it1(Lambda.begin()), itend(Lambda.end());
	   it1 != itend; ++it1, ++row)
	{
	  std::list<size_type> indices;
	  std::list<double> entries;
	  
	  size_type column = 0;
	  for (typename std::set<Index>::const_iterator it2(Lambda.begin());
	     it2 != itend; ++it2, ++column)
	    {
	      double entry = integrate(basis, *it2, *it1);

	      if (entry != 0) {
		indices.push_back(column);
		entries.push_back(entry);
	      }
	    }
	  A_Lambda.set_row(row, indices, entries);
	} 

      // solve A_Lambda*x = b
      Vector<double> b(Lambda.size());
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	b[row] = coeffs.get_coefficient(*it);

      Vector<double> x(b);
      unsigned int iterations;
      CG(A_Lambda, b, x, 1e-12, 200, iterations);
  
      coeffs.clear();
      row = 0;
      for (typename std::set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());
	   it != itend; ++it, ++row)
	coeffs.set_coefficient(*it, x[row]);
#endif
    }
  }

}

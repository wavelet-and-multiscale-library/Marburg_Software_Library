// implementation for cached_problem.h

#include <cmath>
#include <algebra/vector.h>
#include <numerics/eigenvalues.h>

namespace WaveletTL
{
  template <class PROBLEM>
  CachedProblem<PROBLEM>::CachedProblem(const PROBLEM& P)
    : PROBLEM(P)
  {
  }

  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::a(const typename WaveletBasis::Index& lambda,
			    const typename WaveletBasis::Index& nu) const
  {
    // first check whether the lambda-th column already exists in the cache
    typename ColumnCache::iterator col_lb(entries_cache.lower_bound(lambda));
    typename ColumnCache::iterator col_it(col_lb);
    if (col_lb == entries_cache.end() ||
	entries_cache.key_comp()(lambda, col_lb->first))
      {
	// insert a new cache column
	typedef typename ColumnCache::value_type value_type;
	col_it = entries_cache.insert(col_lb, value_type(lambda, Column()));
      }
    
    // for simplicity, extract the column
    Column& col(col_it->second);

    double r = 0;

    // check whether the entry has already been calculated
    typename Column::iterator lb(col.lower_bound(nu));
    typename Column::iterator it(lb);
    if (lb == col.end() ||
	col.key_comp()(nu, lb->first))
      {
	// compute the entry ...
	r = PROBLEM::a(lambda, nu);
	// ... and insert it in the cache
	typedef typename Column::value_type value_type;
	it = col.insert(lb, value_type(nu, r));
      }
    else
      r = it->second;
    
    return r;
  }

  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::norm_A() const
  {
    static double normA = 0.0;
    
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
      
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
    }

    return normA;
  }
   
  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::norm_Ainv() const
  {
    static double normAinv = 0.0;
    
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
      
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
    }

    return normAinv;
  }
}

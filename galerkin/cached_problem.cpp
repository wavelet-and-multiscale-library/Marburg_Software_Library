// implementation for cached_problem.h

#include <cmath>
#include <algebra/vector.h>
#include <numerics/eigenvalues.h>

namespace WaveletTL
{
  template <class PROBLEM>
  CachedProblem<PROBLEM>::CachedProblem(const PROBLEM* P,
					const double estnormA,
					const double estnormAinv)
    : problem(P), normA(estnormA), normAinv(estnormAinv)
  {
  }

  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::a(const Index& lambda,
			    const Index& nu) const
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
	r = problem->a(lambda, nu);
	// ... and insert it into the cache
	typedef typename Column::value_type value_type;
	it = col.insert(lb, value_type(nu, r));
      }
    else
      r = it->second;
    
    return r;
  }
  
  template <class PROBLEM>
  void
  CachedProblem<PROBLEM>::add_level (const Index& lambda,
				     InfiniteVector<double, Index>& w, const int j,
				     const double factor,
				     const int J,
				     const CompressionStrategy strategy) const
  {
    if (problem->local_operator()) {
      typedef std::list<Index> IntersectingList;

      typename StructureCache::iterator col_lb(stiffStructure.lower_bound(lambda));
      typename StructureCache::iterator col_it(col_lb);

      if (col_lb == stiffStructure.end() ||
	  stiffStructure.key_comp()(lambda,col_lb->first))
	{
	  // insert a new column
	  typedef typename StructureCache::value_type value_type;
	  col_it = stiffStructure.insert(col_lb, value_type(lambda, IndexCache()));
	}

      IndexCache& col(col_it->second);
      
      // check wether the level has already been calculated
      typename IndexCache::iterator lb(col.lower_bound(j));
      typename IndexCache::iterator it(lb);
      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
	  IntersectingList nus;
	  // structure an indices for this level have to be
	  // computed
	  intersecting_wavelets(basis(), lambda,
				std::max(j, basis().j0()),
				j == (basis().j0()-1),
				nus);
	  
	  // add the new Indices to our StructureCache
	  typedef typename IndexCache::value_type value_type;
	  it = col.insert(lb, value_type(j, nus));

	  // do the rest of the job
	  const double d1 = problem->D(lambda);
	  if (strategy == St04a) {
	    for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
		 it != itend; ++it) {
	      if (abs(lambda.j()-j) <= J/((double) problem->space_dimension) ||
		  intersect_singular_support(problem->basis(), lambda, *it)) {
		const double entry = a(*it, lambda) / (d1*problem->D(*it));
		w.add_coefficient(*it, entry * factor);
	      }
	    }
	  }
	  else if (strategy == CDD1) {
	    for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
		 it != itend; ++it) {
	      const double entry = a(*it, lambda) / (d1*problem->D(*it));
	      w.add_coefficient(*it, entry * factor);
	    }
	  }

	}
      else {
	// level already exists --> extract level from cache
	const IntersectingList& nus = it->second;

	// do the rest of the job
	const double d1 = problem->D(lambda);
	if (strategy == St04a) {
	  for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	       it != itend; ++it) {
	    if (abs(lambda.j()-j) <= J/((double) problem->space_dimension) ||
		intersect_singular_support(problem->basis(), lambda, *it)) {
	      const double entry = a(*it, lambda) / (d1*problem->D(*it));
	      w.add_coefficient(*it, entry * factor);
	    }
	  }
	}
	else if (strategy == CDD1) {
	  for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	       it != itend; ++it) {
	    const double entry = a(*it, lambda) / (d1*problem->D(*it));
	    w.add_coefficient(*it, entry * factor);
	  }
	  
	}
      }// end else
    }// end if local operator
    else {
      // TODO
    }
  }


  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::norm_A() const
  {
    if (normA == 0.0) {
      cout << "CachedProblem()::norm_A() called..." << endl;

      std::set<Index> Lambda;
      const int j0 = problem->basis().j0();
      const int jmax = j0+1;
      for (Index lambda = problem->basis().first_generator(j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == problem->basis().last_wavelet(jmax)) break;
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

      cout << "... done!" << endl;
    }

    return normA;
  }
   
  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
      cout << "CachedProblem()::norm_Ainv() called..." << endl;

      std::set<Index> Lambda;
      const int j0 = problem->basis().j0();
      const int jmax = j0+1;
      for (Index lambda = problem->basis().first_generator(j0);; ++lambda) {
	Lambda.insert(lambda);
	if (lambda == problem->basis().last_wavelet(jmax)) break;
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

      cout << "... done!" << endl;
    }

    return normAinv;
  }
}

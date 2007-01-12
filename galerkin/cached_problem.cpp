// implementation for cached_problem.h

#include <cmath>
#include <algebra/vector.h>
#include <numerics/eigenvalues.h>
#include <galerkin/galerkin_utils.h>

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
    double r = 0;
    
    if (problem->local_operator()) {

      // BE CAREFUL: KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
      typedef typename Index::type_type generator_type;
      int j = (lambda.e() == generator_type()) ? (lambda.j()-1) : lambda.j();
      
      // check wether entry has already been computed
      typedef std::list<Index> IntersectingList;

      // search for column 'mu'
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu));
      typename ColumnCache::iterator col_it(col_lb);
      
      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(nu, col_lb->first))
	
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(nu, Column()));
	}
	  
      Column& col(col_it->second);
      
      // check wether the level 'lambda' belongs to has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);
      
      

      // no entries have ever been computed for this column and this level
      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
	  // compute whole level block
	  // #### ONLY CDD COMPRESSION STRATEGY IMPLEMENTED ####
	  // #### MAYBE WE ADD TRUNK FOR STEVENSON APPROACH ####
	  // #### LATER.                                    ####
	  
	  // insert a new level
	  typedef typename Column::value_type value_type;
	  it = col.insert(lb, value_type(j, Block()));

	  Block& block(it->second);	  


	  IntersectingList nus;
	   
	  intersecting_wavelets(basis(), nu,
				std::max(j, basis().j0()),
				j == (basis().j0()-1),
				nus);
	  // compute entries
	  for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	       it != itend; ++it) {
	    const double entry = problem->a(*it, nu);
	    typedef typename Block::value_type value_type_block;
	    if (entry != 0.) {
	      block.insert(block.end(), value_type_block(*it, entry));
	      if (Index(*it) == lambda) {
		r = entry;
	      }
	    }
	  } 
	}
      // level already exists --> extract row corresponding to 'lambda'
      else {
	Block& block(it->second);

 	typename Block::iterator block_lb(block.lower_bound(lambda));
 	typename Block::iterator block_it(block_lb);
	// level exists, but in row 'lambda' no entry is available ==> entry must be zero
	if (block_lb == block.end() ||
	    block.key_comp()(lambda, block_lb->first))
	  {
	    r = 0;
	  }
	else {
	  r = block_it->second;
	}
      }
    }
    else {// TODO
      
    }

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
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(lambda));
      typename ColumnCache::iterator col_it(col_lb);
	 
      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(lambda,col_lb->first))

	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(lambda, Column()));
	}
	  
      Column& col(col_it->second);
	  
      // check wether the level has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);

      // no entries have ever been computed for this column and this level
      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{

	  // insert a new level
	  typedef typename Column::value_type value_type;
	  it = col.insert(lb, value_type(j, Block()));

	  IntersectingList nus;
	 
	  intersecting_wavelets(basis(), lambda,
				std::max(j, basis().j0()),
				j == (basis().j0()-1),
				nus);

	  Block& block(it->second);

	  // do the rest of the job
	  const double d1 = problem->D(lambda);
	  if (strategy == St04a) {
	    for (typename IntersectingList::iterator it(nus.begin()), itend(nus.end());
		 it != itend; ++it) {
	      if (abs(lambda.j()-j) <= J/((double) problem->space_dimension) ||
		  intersect_singular_support(problem->basis(), lambda, *it)) {
		const double entry = problem->a(*it, lambda);
		typedef typename Block::value_type value_type_block;

		block.insert(block.end(), value_type_block(*it, entry));
		w.add_coefficient(*it, (entry / (d1*problem->D(*it))) * factor);
	      }
	    }
	  }
	  else if (strategy == CDD1) {
	    for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
		 it != itend; ++it) {
	      const double entry = problem->a(*it, lambda);
	      typedef typename Block::value_type value_type_block;
	      if (entry != 0.)
		block.insert(block.end(), value_type_block(*it, entry));
	      w.add_coefficient(*it, (entry / (d1 * problem->D(*it))) * factor);
	    }
	  }
	}
      else { 
	// level already exists --> extract level from cache
	Block& block(it->second);
	    
	const double d1 = problem->D(lambda);
	// do the rest of the job
	if (strategy == St04a) {
	  for (typename Block::const_iterator it(block.begin()), itend(block.end());
	       it != itend; ++it) {
	    if (abs(lambda.j()-j) <= J/((double) problem->space_dimension) ||
		intersect_singular_support(problem->basis(), lambda, it->first)) {
	      w.add_coefficient(it->first, (it->second / (d1*problem->D(it->first))) * factor);
	    }
	  }
	}
	else if (strategy == CDD1) {
	  for (typename Block::const_iterator it(block.begin()), itend(block.end());
	       it != itend; ++it) {
	    w.add_coefficient(it->first, (it->second / (d1 * problem->D(it->first)))  * factor);
	  }
	}
      }// end else
    }
    else { // TODO
	  
    }
  }


  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::norm_A() const
  {
    if (normA == 0.0) {
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "CachedProblem()::norm_A() called..." << endl;
#endif

      std::set<Index> Lambda;
      const int j0 = problem->basis().j0();
      const int jmax = j0+2;
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

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "... done!" << endl;
#endif
    }

    return normA;
  }
   
  template <class PROBLEM>
  double
  CachedProblem<PROBLEM>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "CachedProblem()::norm_Ainv() called..." << endl;
#endif

      std::set<Index> Lambda;
      const int j0 = problem->basis().j0();
      const int jmax = j0+2;
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

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "... done!" << endl;
#endif
    }

    return normAinv;
  }


  template <class PROBLEM>
  CachedProblemFromFile<PROBLEM>::CachedProblemFromFile
  (const PROBLEM* P,
   const char* filename,
   const int jmax,
   const double estnormA,
   const double estnormAinv)
    : problem(P), jmax_(jmax), normA(estnormA), normAinv(estnormAinv)
  {
    entries_cache.matlab_input(filename);
  }

  template <class PROBLEM>
  inline
  double
  CachedProblemFromFile<PROBLEM>::a(const Index& lambda,
				    const Index& nu) const
  {
    return entries_cache.get_entry(lambda.number(), nu.number());
  }

  template <class PROBLEM>
  void
  CachedProblemFromFile<PROBLEM>::add_level (const Index& lambda,
					     InfiniteVector<double, Index>& w, const int j,
					     const double factor,
					     const int J,
					     const CompressionStrategy strategy) const
  {
    if (problem->local_operator()) {
      typedef std::list<Index> IntersectingList;
      IntersectingList nus;
      intersecting_wavelets(basis(), lambda,
			    std::max(j, basis().j0()),
			    j == (basis().j0()-1),
			    nus);

      // do the rest of the job, we only implement the CDD1 compression strategy here!!!
      const double d1 = problem->D(lambda);
      for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	   it != itend; ++it) {
	const double entry = problem->a(*it, lambda);
	w.add_coefficient(*it, (entry / (d1 * problem->D(*it))) * factor);
      }
    } else {
      // integral operator case, TODO
    }
  }

  template <class PROBLEM>
  double
  CachedProblemFromFile<PROBLEM>::norm_A() const
  {
    if (normA == 0.0) {
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "CachedProblemFromFile()::norm_A() called..." << endl;
#endif

      std::set<Index> Lambda;
      const int j0 = problem->basis().j0();
      const int jmax = j0+2;
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

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "... done!" << endl;
#endif
    }

    return normA;
  }
   
  template <class PROBLEM>
  double
  CachedProblemFromFile<PROBLEM>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "CachedProblemFromFile()::norm_Ainv() called..." << endl;
#endif

      std::set<Index> Lambda;
      const int j0 = problem->basis().j0();
      const int jmax = j0+2;
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

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "... done!" << endl;
#endif
    }

    return normAinv;
  }

}

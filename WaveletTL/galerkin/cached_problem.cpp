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
        
        
      const int lambda_num = lambda.number();
      //cout << "number: " << lambda_num << endl;
      const int nu_num = nu.number();
      //cout << "number2: " << nu_num << endl;
      
      // BE CAREFUL: KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
      typedef typename Index::type_type generator_type;
      int j = (lambda.e() == generator_type()) ? (lambda.j()-1) : lambda.j();
      
      // check wether entry has already been computed
      typedef std::list<Index> IntersectingList;

      // search for column 'mu'
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
      typename ColumnCache::iterator col_it(col_lb);
      
      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(nu_num, col_lb->first))
	
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
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
            //cout << *it << ", " << nu << ": " << entry << endl;
	    typedef typename Block::value_type value_type_block;
	    if (entry != 0.) {
	      block.insert(block.end(), value_type_block((*it).number(), entry));
	      if ((*it).number() == lambda_num) {
		r = entry;
	      }
	    }
	  } 
	}
      // level already exists --> extract row corresponding to 'lambda'
      else {
	Block& block(it->second);

 	//typename Block::iterator block_lb(block.lower_bound(lambda));
	typename Block::iterator block_lb(block.lower_bound(lambda_num));
 	typename Block::iterator block_it(block_lb);
	// level exists, but in row 'lambda' no entry is available ==> entry must be zero
	if (block_lb == block.end() ||
	    block.key_comp()(lambda_num, block_lb->first))
	  {
	    r = 0;
	  }
	else {
	  r = block_it->second;
	}
      }
    }
    else {
      // for nonlocal operators, we put full level blocks into the cache, regardless of support intersections

      const int lambda_num = lambda.number();
      const int nu_num = nu.number();
      
      // BE CAREFUL: KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
      typedef typename Index::type_type generator_type;
      int j = (lambda.e() == generator_type()) ? (lambda.j()-1) : lambda.j();
      
      // search for column 'mu'
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
      typename ColumnCache::iterator col_it(col_lb);
      
      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(nu_num, col_lb->first))
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
	}
      
      Column& col(col_it->second);
      
      // check wether the level block which 'lambda' belongs to has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);
      
      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
	  // no entries have ever been computed for this column and this level
	  
	  // compute whole level block
	  
	  // insert a new level
	  typedef typename Column::value_type value_type;
	  it = col.insert(lb, value_type(j, Block()));

	  Block& block(it->second);	  

	  // collect all indices in the level block
	  typedef std::list<Index> IndexList;
	  IndexList nus;
	  if (j == (basis().j0()-1)) {
	    // generators on level j0
	    for (Index lambda1 = basis().first_generator(basis().j0());; ++lambda1) {
	      nus.push_back(lambda1);
	      if (lambda1 == basis().last_generator(basis().j0())) break;
	    }
	  } else {
	    // wavelets on level j
	    for (Index lambda1 = basis().first_wavelet(j);; ++lambda1) {
	      nus.push_back(lambda1);
	      if (lambda1 == basis().last_wavelet(j)) break;
	    }
	  }

	  // compute entries
	  for (typename IndexList::const_iterator it(nus.begin()), itend(nus.end());
	       it != itend; ++it) {
	    const double entry = problem->a(*it, nu);
	    typedef typename Block::value_type value_type_block;
	    if (fabs(entry) > 1e-16 ) {
	      block.insert(block.end(), value_type_block((*it).number(), entry));
	      if ((*it).number() == lambda_num) {
		r = entry;
	      }
	    }
	  } 
	}
      // level already exists --> extract row corresponding to 'lambda'
      else {
	Block& block(it->second);
	
	typename Block::iterator block_lb(block.lower_bound(lambda_num));
 	typename Block::iterator block_it(block_lb);
	// level exists, but in row 'lambda' no entry is available ==> entry must be zero
	if (block_lb == block.end() ||
	    block.key_comp()(lambda_num, block_lb->first)) {
	  r = 0;
	}
	else {
	  r = block_it->second;
	}
      }
    }
    
    return r;
  }
  
  template <class PROBLEM>
  void
  CachedProblem<PROBLEM>::add_level(const Index& lambda,
				     //InfiniteVector<double, Index>& w,
				     Vector<double>& w,
				     const int j,
				     const double factor,
				     const int J,
				     const CompressionStrategy strategy,
                                     const int jmax,
                                     const int pmax,
                                     const double a,
                                     const double b) const
  {
    if (problem->local_operator()) {

      int lambda_num = lambda.number();

      // search for column 'lambda'
      typedef std::list<Index> IntersectingList;
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(lambda_num));
      typename ColumnCache::iterator col_it(col_lb);
	 
      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(lambda_num,col_lb->first))
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(lambda_num, Column()));
	}

      Column& col(col_it->second);

      // check wether the level has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);

      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
	  // no entries have ever been computed for this column and this level
	  
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
	    for (typename IntersectingList::iterator it2(nus.begin()), itend2(nus.end());
		 it2 != itend2; ++it2) {
	      if (abs(lambda.j()-j) <= J/((double) problem->space_dimension) ||
		  intersect_singular_support(problem->basis(), lambda, *it2)) {
		const double entry = problem->a(*it2, lambda);
		typedef typename Block::value_type value_type_block;
		block.insert(block.end(), value_type_block((*it2).number(), entry));
		//w.add_coefficient(*it2, (entry / (d1*problem->D(*it2))) * factor);
		w[(*it2).number()] += (entry / (d1*problem->D(*it2))) * factor;
	      }
	    }
	  }
	  else if (strategy == CDD1 || strategy == DKOR) {
	    for (typename IntersectingList::const_iterator it2(nus.begin()), itend2(nus.end());
		 it2 != itend2; ++it2) {
	      const double entry = problem->a(*it2, lambda);
	      typedef typename Block::value_type value_type_block;
	      if (entry != 0.) {
		block.insert(block.end(), value_type_block((*it2).number(), entry));
	      //w.add_coefficient(*it2, (entry / (d1 * problem->D(*it2))) * factor);
		w[(*it2).number()] += (entry / (d1*problem->D(*it2))) * factor;
	      }
	    }
	  }
	}
      else {
	// level already exists --> extract level from cache

	Block& block(it->second);
	
	const double d1 = problem->D(lambda);
	
	// do the rest of the job
	if (strategy == St04a) {
	  
	  for (typename Block::const_iterator it2(block.begin()), itend2(block.end());
	       it2 != itend2; ++it2) {
 	    if (abs(lambda.j()-j) <= J/((double) problem->space_dimension) ||
 		intersect_singular_support(problem->basis(), lambda, *(problem->basis().get_wavelet(it2->first)))) {
// 	      w.add_coefficient(*(problem->basis().get_wavelet(it2->first)),
// 				(it2->second / (d1*problem->D(*(problem->basis().get_wavelet(it2->first))))) * factor);
	      
	      w[it2->first] += (it2->second / (d1*problem->D(*(problem->basis().get_wavelet(it2->first))))) * factor;

	    }
	  }
	}
	else if (strategy == CDD1 || strategy == DKOR) {
	  for (typename Block::const_iterator it2(block.begin()), itend2(block.end());
	       it2 != itend2; ++it2) {
// 	    w.add_coefficient(*(problem->basis().get_wavelet(it2->first)),
// 			      (it2->second / (d1 * problem->D( *(problem->basis().get_wavelet(it2->first)) )))  * factor);
	    w[it2->first] += (it2->second / (d1*problem->D(*(problem->basis().get_wavelet(it2->first))))) * factor;
	  }
	}
      }// end else
    }
    else {
      // for nonlocal operators, we put full level blocks into the cache, regardless of support intersections

      const int lambda_num = lambda.number();

      // search for column 'lambda'
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(lambda_num));
      typename ColumnCache::iterator col_it(col_lb);
      
      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(lambda_num, col_lb->first))
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(lambda_num, Column()));
	}
      
      Column& col(col_it->second);
      
      // check wether the level block which 'lambda' belongs to has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);
      
      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
	  // no entries have ever been computed for this column and this level
	  
	  // compute whole level block
	  
	  // insert a new level
 	  typedef typename Column::value_type value_type;
 	  it = col.insert(lb, value_type(j, Block()));

 	  Block& block(it->second);  

	  // collect all indices in the level block
	  typedef std::list<Index> IndexList;
	  IndexList nus;
	  if (j == (basis().j0()-1)) {
	    // generators on level j0
	    for (Index lambda1 = basis().first_generator(basis().j0());; ++lambda1) {
	      nus.push_back(lambda1);
	      if (lambda1 == basis().last_generator(basis().j0())) break;
	    }
	  } else {
	    // wavelets on level j
	    for (Index lambda1 = basis().first_wavelet(j);; ++lambda1) {
	      nus.push_back(lambda1);
	      if (lambda1 == basis().last_wavelet(j)) break;
	    }
	  }

	  // compute entries
	  const double d1 = problem->D(lambda);
	  for (typename IndexList::const_iterator it2(nus.begin()), itend2(nus.end());
	       it2 != itend2; ++it2) {
	    const double entry = problem->a(*it2, lambda);
 	    typedef typename Block::value_type value_type_block;
	    if (fabs(entry) > 1e-16) {
  	      block.insert(block.end(), value_type_block((*it2).number(), entry));
	      w[(*it2).number()] += (entry / (d1*problem->D(*it2))) * factor;
	    }
	  } 
	}
      else {
 	// level already exists --> extract level from cache

 	Block& block(it->second);

 	const double d1 = problem->D(lambda);

 	for (typename Block::const_iterator it(block.begin()), itend(block.end());
 	     it != itend; ++it) {
 	  w[it->first] += (it->second / (d1*problem->D(*(problem->basis().get_wavelet(it->first))))) * factor;
 	}
      }
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
        //if (i==7) break;
      }
      SparseMatrix<double> A_Lambda;
      
      Matrix<double> evecs;
      Vector<double> evals;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      A_Lambda.compress(1e-10);
      
#if 0
      double help;
      unsigned int iterations;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      normAinv = 1./help;
      cout << "Norm:" << normA << endl;
      cout << "inv Norm:" << normAinv << endl;
      //cout << help << endl;
      cout << "Iter: " << iterations << endl;
#else
      //LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      SymmEigenvalues(A_Lambda, evals, evecs);
      
      cout << "Eigenwerte: " << evals << endl;
      //cout << "Eigenvektoren: " << endl << evecs << endl;
      normA = evals(evals.size()-1);
      cout << "normA: " << normA << endl;
      int i = 0;
      while(abs(evals(i))<1e-10){
          ++i;
      }
      cout << "rezinvers: " <<evals(i) << endl;;
      normAinv = 1./evals(i);
      cout << "inversnorm: " << normAinv << endl;
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
        //cout << lambda << endl;
	if (lambda == problem->basis().last_wavelet(jmax)) break;
      }
      //cout << "Schritt 1" << endl;
      SparseMatrix<double> A_Lambda;
      Matrix<double> evecs;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      //cout << "Matrix aufgestellt" << endl;
      //cout << A_Lambda << endl;
      Vector<double> evals;
      
#if 1
      cout << "Beginne Normberechnung" << endl;
      SymmEigenvalues(A_Lambda, evals, evecs);
      normA = evals(evals.size()-1);
      int i = 0;
      while(abs(evals(i))<1e-10){
          ++i;
      }
      cout << "rezinvers: " <<evals(i) << endl;;
      normAinv = 1./evals(i);
      cout << "Norm " << normA << endl;
      cout << "Inv Norm " << normAinv << endl;
      
#else
      Vector<double> xk(Lambda.size(), false);
      xk = 1;
      unsigned int iterations;
      normAinv = 1./InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
      cout << "Inv Norm " << normAinv << endl;
      cout << iterations << endl;
      //cout << A_Lambda << endl;
      
      //cout << "Eigenvektoren: " << evecs << endl;
#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "... done!" << endl;
#endif
    }

    return normAinv;
  }

  // ############## THE FOLLOWING TWO ROUTINES ARE PURELY EXPERIMENTAL AT THE MOMENT! ###########

//     // determining the first index of a new level in the index set
//     list<int> level_pos;
//     typename std::set<int>::const_iterator win_it = window.begin();
//     level_pos.push_back(*win_it);

//     while ( true ) {
//       const Index* ind = problem->basis().get_wavelet(*win_it);
//       int old_lev = (( ind->e() == generator_type()) ? ( ind->j()-1) : ind->j());
//       win_it++;
//       if ( win_it == window.end() )
// 	break;
//       const Index* ind2 = problem->basis().get_wavelet(*win_it);
//       int new_lev = (( ind2->e() == generator_type()) ? ( ind2->j()-1) : ind2->j());
//       if (new_lev != old_lev)
// 	level_pos.push_back(*win_it);
//     }
    
// 	typename list<int>::const_iterator level_pos_it = level_pos.begin();
// 	level_pos_it++; // starting position of second level occuring

  template <class PROBLEM>
  void
  CachedProblem<PROBLEM>::apply(const std::set<int>& window, const Vector<double>& x,
				Vector<double>& res) const
  {
//     for (typename std::set<int>::const_iterator iter = window.begin(); iter != window.end(); iter++) {
//       cout << *(problem->basis().get_wavelet(*iter)) << endl;
//     }
//     cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    res.resize(x.size());
    typedef typename Index::type_type generator_type;
    // cout << " size = " << entries_cache.size() << endl;  
    if (entries_cache.size() == 0) {
      unsigned int l = 0;
      for (typename std::set<int>::const_iterator win_it_col = window.begin(); win_it_col != window.end(); win_it_col++, l++) {
	const double d1 = problem->D(*(problem->basis().get_wavelet(*win_it_col)));
	unsigned int k = 0;
	for (typename std::set<int>::const_iterator win_it_row = window.begin(); win_it_row != window.end(); win_it_row++, k++) {
	  res[k] += x[l] * this->a(*(problem->basis().get_wavelet(*win_it_row)),*(problem->basis().get_wavelet(*win_it_col)))
	    / (d1*problem->D(*(problem->basis().get_wavelet(*win_it_row))));
	}
      }
      return;
    }

    typename std::set<int>::const_iterator win_it_col = window.begin();
    unsigned int l = 0;

    for (typename ColumnCache::iterator col_it = entries_cache.begin();
	 col_it != entries_cache.end(); col_it++) {
      if (col_it->first < *win_it_col) {
	continue;
      }
      else if (col_it->first > *win_it_col){
	//	cout << "inserting column " << *win_it_col << endl;
	const double d1 = problem->D(*(problem->basis().get_wavelet(*win_it_col)));
	unsigned int k = 0;
	for (typename std::set<int>::const_iterator win_it_row = window.begin();
	     win_it_row != window.end(); win_it_row++, k++) {
	  res[k] += x[l] * this->a(*(problem->basis().get_wavelet(*win_it_row)),
			     *(problem->basis().get_wavelet(*win_it_col)))
	    / (d1*problem->D(*(problem->basis().get_wavelet(*win_it_row))));
	}
	col_it--;
// 	cout << "column should be " << (*win_it_col) << " and it is " << col_it->first << endl;
// 	cout << "done inserting column " << (*win_it_col) << endl;
	win_it_col++;
	l++;
	if(l == res.size())
	  return;
      }
      else if (col_it->first == *win_it_col) {
	typename std::set<int>::const_iterator win_it_row = window.begin();
	unsigned int k = 0;
	
	for (typename Column::const_iterator block_it = (col_it->second).begin();
	     block_it != (col_it->second).end(); block_it++) {

	  const Index* rowind = problem->basis().get_wavelet(*win_it_row);

	  int j = (rowind->e() == generator_type()) ? (rowind->j()-1) : rowind->j();

	  if (block_it->first < j) {
	    continue;
	  }
	  else if (block_it->first > j) { // insert this block!
	    //	    cout << "inserting level " << j << endl;
	    Column& col(col_it->second);
	    typename Column::iterator lb(col.lower_bound(j));
	    typename Column::iterator it(lb);
	     
	    typedef typename Column::value_type value_type;
	    it = col.insert(lb, value_type(j, Block()));
	    std::list<Index> nus;

	    intersecting_wavelets(basis(), *(problem->basis().get_wavelet(col_it->first)),
				  std::max(j, basis().j0()),
				  j == (basis().j0()-1),
				  nus);
	    Block& block(it->second);
	    const double d1 = problem->D(*(problem->basis().get_wavelet(col_it->first)));
	     
	    for (typename std::list<Index>::const_iterator it(nus.begin()), itend(nus.end());
		 it != itend; ++it) {
	      //	      cout << "nonzero entry at " << (*it).number() << endl;
	      while (*win_it_row < (*it).number() && win_it_row != window.end()) {
		win_it_row++;
		k++;
	      }
	      //	      cout << "k=" << k << endl;
	      const double entry = problem->a(*it, problem->basis().get_wavelet(col_it->first));
	      typedef typename Block::value_type value_type_block;
	      if (entry != 0.)
		block.insert(block.end(), value_type_block((*it).number(), entry));
	       
	      // at this point holds either (*win_it_row).number() < (*it).number() or (*win_it_row).number() == (*it).number()
	      if ( *win_it_row == (*it).number() ) {
		if (entry != 0.)
		  res[k] += x[l] * (entry / (d1*problem->D(*it)));
	      }
	    }
	    if (win_it_row == window.end())
	      break;
	    block_it--;
	    // leaving this level
	    while ( (( problem->basis().get_wavelet(*win_it_row)->e() == generator_type()) ?
		     ( problem->basis().get_wavelet(*win_it_row)->j()-1) :
		     problem->basis().get_wavelet(*win_it_row)->j()) == j
		    && win_it_row != window.end()) {
	      win_it_row++;
	      k++;
	    }
	    if (win_it_row == window.end())
	      break;
	  }
	  else if (block_it->first == j) { // extract existing block!
	    //	    cout << "extracting level " << j << " in column " << col_it->first << endl;
	    const double d1 = problem->D(*(problem->basis().get_wavelet(col_it->first)));
	    for (typename Block::const_iterator it = block_it->second.begin(); it != block_it->second.end(); ++it) {
	      //	      cout << "nonzero entry in row " << it->first << endl;
	      const Index* ind = problem->basis().get_wavelet(it->first);
	      while (*win_it_row < ind->number() && k < res.size()) {
		win_it_row++;
		k++;
	      }
	      //	      cout << "k = " << k << endl;
	      if (k == res.size())
		break;
	      // at this point holds either (*win_it_row).number() < (*it).number() or (*win_it_row).number() == (*it).number()
	      if ( *win_it_row == ind->number() ) {
		//		cout << "putting row " << *win_it_row << endl;
		res[k] += x[l] * (it->second / (d1*problem->D(*ind)));
	      }
	    
	    }// end loop over level block
	    if (k == res.size())
	      break;
	    //cout << "leaving the level " << endl;
	    // leaving this level
	    while ( (( problem->basis().get_wavelet(*win_it_row)->e() == generator_type()) ?
		     ( problem->basis().get_wavelet(*win_it_row)->j()-1) :
		     problem->basis().get_wavelet(*win_it_row)->j()) == j
		    && k < res.size()) {
	      //	      cout << "kk = " << k << endl;
	      win_it_row++;
	      k++;
	    }
	    if (k == res.size())
	      break;
	  }// end else if (block_it->first == j)
	}// end loop over level blocks
	win_it_col++;
	l++;
	if (l == res.size())
	  return;
      }// end if existing column found
    }// end loop over all columns in cache

    // insert the remaining columns
    while (l < res.size()) {
      //      cout << "inserting column " << *win_it_col << endl;
      const double d1 = problem->D(*(problem->basis().get_wavelet(*win_it_col)));
      unsigned int k = 0;
      for (typename std::set<int>::const_iterator win_it_row = window.begin();
	   win_it_row != window.end(); win_it_row++, k++) {
	res[k] += x[l] * this->a(*(problem->basis().get_wavelet(*win_it_row)),
				 *(problem->basis().get_wavelet(*win_it_col)))
	  / (d1*problem->D(*(problem->basis().get_wavelet(*win_it_row))));
      }
      win_it_col++;
      l++;
    }
  }
  
  template <class PROBLEM>
  bool
  CachedProblem<PROBLEM>::CG(const std::set<int>& window, const Vector<double> &b, Vector<double> &xk,
	  const double tol, const unsigned int maxiter, unsigned int& iterations)
  {
    // see: "Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods"

    Vector<double> rk(window.size(), false),
      zk(window.size(), false),
      pk(window.size(), false),
      Apk(window.size(), false);

    // first (negative) residual
    //cout << "entering windowed matrix vector product..." << endl;
    apply(window, xk, rk);
    //cout << "done windowed matrix vector product..." << endl;
    //cout << rk << endl;
    rk.subtract(b);
    const double normr0 = l2_norm_sqr(rk);
    double normrk = normr0, rhok = 0, oldrhok = 0;
    for (iterations = 1; normrk/normr0 > tol*tol && iterations <= maxiter; iterations++)
      {
	//P.apply_preconditioner(rk, zk);
	zk = rk;
	rhok = rk * zk;

	if (iterations == 1) // TODO: shift this case upwards!
	  pk = zk;
	else
	  pk.sadd(rhok/oldrhok, zk);
	//cout << "entering windowed matrix vector product..." << endl;
	apply(window, pk, Apk);
	
	//cout << "done windowed matrix vector product..." << endl;
	//cout << Apk << endl;
	const double alpha = rhok/(pk*Apk);
	xk.add(-alpha,  pk);
	rk.add(-alpha, Apk);
	normrk = l2_norm_sqr(rk);

	oldrhok = rhok;
      }

    return (iterations <= maxiter);
  }


  //
  //
  // CachedProblemLocal

  template <class PROBLEM>
  CachedProblemLocal<PROBLEM>::CachedProblemLocal(const PROBLEM* P,
						  const double estnormA,
						  const double estnormAinv)
    : problem(P), normA(estnormA), normAinv(estnormAinv)
  {
    entries_cache.resize(basis().n_p());
  }

  template <class PROBLEM>
  double
  CachedProblemLocal<PROBLEM>::a(const Index& lambda,
				 const Index& nu) const
  {
    const int p = lambda.p();
     
    double r = 0;
    
    if (problem->local_operator()) {
      const int lambda_num = lambda.number();
      const int nu_num = nu.number();


      // BE CAREFUL: KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
      typedef typename Index::type_type generator_type;
      int j = (lambda.e() == generator_type()) ? (lambda.j()-1) : lambda.j();
      
      // check wether entry has already been computed
      typedef std::list<Index> IntersectingList;

      // search for column 'mu'
      typename ColumnCache::iterator col_lb(entries_cache[p].lower_bound(nu_num));
      typename ColumnCache::iterator col_it(col_lb);
      
      if (col_lb == entries_cache[p].end() ||
	  entries_cache[p].key_comp()(nu_num, col_lb->first))
	
	{
#ifdef ONE_D
	  return problem->a(lambda, nu);
#else
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache[p].insert(col_lb, value_type(nu_num, Column()));
#endif
	}
	  
      Column& col(col_it->second);
      
      // check wether the level 'lambda' belongs to has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);
      
      

      // no entries have ever been computed for this column and this level
      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
#ifdef ONE_D
	  return problem->a(lambda, nu);
#else
	  // compute whole level block
	  // #### ONLY CDD COMPRESSION STRATEGY IMPLEMENTED ####
	  // #### MAYBE WE ADD TRUNK FOR STEVENSON APPROACH ####
	  // #### LATER.                                    ####
	  
	  // insert a new level
	  typedef typename Column::value_type value_type;
	  it = col.insert(lb, value_type(j, Block()));

	  Block& block(it->second);	  
	  IntersectingList nus;
	  
	  intersecting_wavelets_on_patch(basis(), nu,
					 p,
					 std::max(j, basis().j0()),
					 j == (basis().j0()-1),
					 nus);
	  // compute entries
	  for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
	       it != itend; ++it) {
	    const double entry = problem->a(*it, nu);
	    typedef typename Block::value_type value_type_block;
	    if (entry != 0.) {
	      block.insert(block.end(), value_type_block((*it).number(), entry));
	      if ((*it).number() == lambda_num) {
		r = entry;
	      }
	    }
	  }
#endif	  
	}
      // level already exists --> extract row corresponding to 'lambda'
      else {
	Block& block(it->second);
	
 	//typename Block::iterator block_lb(block.lower_bound(lambda));
	typename Block::iterator block_lb(block.lower_bound(lambda_num));
 	typename Block::iterator block_it(block_lb);
	// level exists, but in row 'lambda' no entry is available ==> entry must be zero
	if (block_lb == block.end() ||
	    block.key_comp()(lambda_num, block_lb->first))

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
  CachedProblemLocal<PROBLEM>::add_level (const Index& lambda, 
					  const int p,
					  Vector<double>& w,
					  const int j,
					  const double factor,
					  const int J,
					  const CompressionStrategy strategy) const
  {
    if (problem->local_operator()) {

      int lambda_num = lambda.number();

      typedef std::list<Index> IntersectingList;
      typename ColumnCache::iterator col_lb(entries_cache[p].lower_bound(lambda_num));
      typename ColumnCache::iterator col_it(col_lb);
	 
      if (col_lb == entries_cache[p].end() ||
	  entries_cache[p].key_comp()(lambda_num,col_lb->first))
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache[p].insert(col_lb, value_type(lambda_num, Column()));
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

	  intersecting_wavelets_on_patch(basis(), lambda,
					 p,
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
		block.insert(block.end(), value_type_block((*it).number(), entry));
		//w.add_coefficient(*it, (entry / (d1*problem->D(*it))) * factor);
		w[(*it).number()] += (entry / (d1*problem->D(*it))) * factor;
	      }
	    }
	  }
	  else if (strategy == CDD1) {
	    for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
		 it != itend; ++it) {
	      const double entry = problem->a(*it, lambda);
	      typedef typename Block::value_type value_type_block;
	      if (entry != 0.)
		block.insert(block.end(), value_type_block((*it).number(), entry));
	      //w.add_coefficient(*it, (entry / (d1 * problem->D(*it))) * factor);
	      w[(*it).number()] += (entry / (d1*problem->D(*it))) * factor;
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
 		intersect_singular_support(problem->basis(), lambda, *(problem->basis().get_wavelet(it->first)))) {
// 	      w.add_coefficient(*(problem->basis().get_wavelet(it->first)),
// 				(it->second / (d1*problem->D(*(problem->basis().get_wavelet(it->first))))) * factor);
	      
	      w[it->first] += (it->second / (d1*problem->D(*(problem->basis().get_wavelet(it->first))))) * factor;

	    }
	  }
	}
	else if (strategy == CDD1) {
	  for (typename Block::const_iterator it(block.begin()), itend(block.end());
	       it != itend; ++it) {
// 	    w.add_coefficient(*(problem->basis().get_wavelet(it->first)),
// 			      (it->second / (d1 * problem->D( *(problem->basis().get_wavelet(it->first)) )))  * factor);
	    w[it->first] += (it->second / (d1*problem->D(*(problem->basis().get_wavelet(it->first))))) * factor;
	  }
	}
      }// end else
    }
    else { // TODO
	  
    }
  }

  template <class PROBLEM>
  double
  CachedProblemLocal<PROBLEM>::norm_A() const
  {
    // TODO
    return 1;;
  }
   
  template <class PROBLEM>
  double
  CachedProblemLocal<PROBLEM>::norm_Ainv() const
  {
    // TODO
    return 1;
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
					     //InfiniteVector<double, Index>& w,
					     Vector<double>& w,
					     const int j,
					     const double factor,
					     const int J,
					     const CompressionStrategy strategy) const
  {
    typedef typename Vector<double>::size_type size_type;

//     if (problem->local_operator()) {
      // remark: the following loop can still be optimized
      Index mu = (j == (basis().j0()-1)
 		  ? basis().first_generator(basis().j0())
 		  : basis().first_wavelet(j));
      size_type start_index = (j == (basis().j0()-1)
 			       ? 0
 			       : basis().Deltasize(j));
      size_type end_index = basis().Deltasize(j+1);
      size_type lambda_number = lambda.number();
//       const double d1 = sqrt(entries_cache.get_entry(lambda_number, lambda_number)); // fails for D==1
      const double d1 = problem->D(lambda);
      for (size_type i = start_index; i < end_index; i++, ++mu) {
	const double entry = entries_cache.get_entry(i,lambda_number);
	if (entry != 0)
// 	  w.add_coefficient(mu, entry * factor /
// // 			    (d1 * sqrt(entries_cache.get_entry(i, i)))); // fails for D==1
// 			    (d1 * problem->D(mu)));
	w[mu.number()] += entry * factor /
// 			    (d1 * sqrt(entries_cache.get_entry(i, i)))); // fails for D==1
	                    (d1 * problem->D(mu));

      }
//     }
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

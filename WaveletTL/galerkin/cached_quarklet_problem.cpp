// implementation for cached_quarklet_problem.h

#include <cmath>
#include <algebra/vector.h>
#include <numerics/eigenvalues.h>
#include <galerkin/galerkin_utils.h>

namespace WaveletTL
{
  template <class PROBLEM>
  CachedQuarkletProblem<PROBLEM>::CachedQuarkletProblem(const PROBLEM* P,
					const double estnormA,
					const double estnormAinv)
    : problem(P), normA(estnormA), normAinv(estnormAinv)
  {
  }
  
  template <class PROBLEM>
  int 
  CachedQuarkletProblem<PROBLEM>::number (const Index& lambda, const int jmax) const{
      
      QuarkletFrame myframe(frame());
      const int quarkletsonzerolevel = (myframe.last_wavelet(jmax)).number()+1;
      const int quarkletsonplevel = ((myframe.last_wavelet(jmax,lambda.p())).number()+1);
      const int lambda_num = (lambda.p()==0) ? lambda.number() : lambda.number() + (lambda.p()-1) * quarkletsonplevel + quarkletsonzerolevel;
      return lambda_num;
  }

  template <class PROBLEM>
  double
  CachedQuarkletProblem<PROBLEM>::a(const Index& lambda,
			    const Index& nu) const
  {
      
    //const int lambda_num = number(lambda,2);
    double r = 0;

    QuarkletFrame myframe(frame());

//!!! number() seems to be a huge bottleneck
    const int jmax = myframe.get_jmax_();
    const int lambda_num = number(lambda, jmax);
//    const int lambda_num = lambda.number();
    const int nu_num = number(nu,jmax);
//    const int nu_num = nu.number();
    
//    const int lambda_num = lambda.number();
//      //cout << "number: " << lambda_num << endl;
//      const int nu_num = nu.number();
//      //cout << "number2: " << nu_num << endl;

    //cout << "Punkt 1" << endl;
    // BE CAREFUL: KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
    typedef typename Index::type_type generator_type;
    int j = (lambda.e() == generator_type()) ? (lambda.j()-1) : lambda.j();
    int p = lambda.p();

    // check wether entry has already been computed
    typedef std::list<Index> IntersectingList;

    // search for column 'nu'
    typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
    typename ColumnCache::iterator col_it(col_lb);
    // no entries have ever been computed for this column
    if (col_lb == entries_cache.end() ||
        entries_cache.key_comp()(nu_num, col_lb->first))

      {
        // insert a new column
        typedef typename ColumnCache::value_type value_type;
        col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
      }
       //cout << "Punkt 2" << endl; 
    Column& col(col_it->second);

    // check wether the polynomial 'lambda' belongs to has already been calculated
    typename Column::iterator block_lb(col.lower_bound(p));
    typename Column::iterator block_it(block_lb);
    // no entries have ever been computed for this column and this polynomial
    if (block_lb == col.end() ||
        col.key_comp()(p, block_lb->first))
      {
        // insert a new polynomial
         typedef typename Column::value_type value_type;
         block_it = col.insert(block_lb, value_type(p, Block()));
      }

    Block& block(block_it->second);

    // check wether the level 'lambda' belongs to has already been calculated
    typename Block::iterator lb(block.lower_bound(j));
    typename Block::iterator it(lb);
    // no entries have ever been computed for this column, polynomial and level
    if (lb == block.end() ||
        block.key_comp()(j, lb->first))
        {
           //cout << "Punkt 3" << endl; 
            // insert a new level
            typedef typename Block::value_type value_type;
            it = block.insert(lb, value_type(j, Subblock()));
      

            Subblock& subblock(it->second);

            IntersectingList nus;


            intersecting_quarklets(frame(), nu,
                                  std::max(j, frame().j0()),
                                  j == (frame().j0()-1),
                                  nus, p);
        //            for (typename IntersectingList::const_iterator it100(nus.begin()), itend(nus.end());
        //                   it100 != itend; ++it100) {
        //                cout << *it100 << endl;
        //            }


            // compute entries
            for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());
                   it != itend; ++it) {
                const double entry = problem->a(*it, nu);
        //                cout << *it << ", " << nu << ": " << entry << endl;
        //                cout << (*it).number()+(*it).p()*quarkletsonplevel << endl;
                typedef typename Subblock::value_type value_type_subblock;
                if (entry != 0.) {
                    subblock.insert(subblock.end(), value_type_subblock(number(*it,jmax), entry));
                    if (number(*it,jmax) == lambda_num) {
                        r = entry;
                    }
                }
            }
          
          //cout << "nu" << nu << endl;
          //cout << block.end() << endl;
	}
      // level already exists --> extract row corresponding to 'lambda'
    else {
        //cout << "ja" << endl;
        Subblock& subblock(it->second);

      //typename Block::iterator block_lb(block.lower_bound(lambda));
      typename Subblock::iterator subblock_lb(subblock.lower_bound(lambda_num));
      typename Subblock::iterator subblock_it(subblock_lb);
      // level exists, but in row 'lambda' no entry is available ==> entry must be zero
      if (subblock_lb == subblock.end() ||
          subblock.key_comp()(lambda_num, subblock_lb->first))
        {
          r = 0;
        }
      else {
        r = subblock_it->second;
      }
    }
    

    return r;
  }
  
  
  
  template <class PROBLEM>
  void
  CachedQuarkletProblem<PROBLEM>::add_level (const Index& lambda,
				     //InfiniteVector<double, Index>& w,
				     Vector<double>& w,
                                     const int p, 
				     const int j,
				     const double factor,
				     const int J,
				     const CompressionStrategy strategy,
                                     const int jmax,
                                     const int pmax,          
                                     const double a,
                                     const double b) const
  {
//      cout << "bin in addlevel drin" << endl;
//      cout << lambda << endl;
//      cout << j << endl;
      
//      QuarkletFrame myframe(frame());
      
      

      //!!! number() seems to be a huge bottleneck
      const int lambda_num = number(lambda, jmax);
//      const int lambda_num = lambda.number();
      //cout << lambda << ": " << lambda_num << endl;
      
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
      
        // check wether the polynomial 'lambda' belongs to has already been calculated
      typename Column::iterator block_lb(col.lower_bound(p));
      typename Column::iterator block_it(block_lb);
      // no entries have ever been computed for this column and this polynomial
      if (block_lb == col.end() ||
          col.key_comp()(p, block_lb->first))
        {
          // insert a new polynomial
           typedef typename Column::value_type value_type;
           block_it = col.insert(block_lb, value_type(p, Block()));
        }

      Block& block(block_it->second);
      
      
        // check wether the level 'lambda' belongs to has already been calculated
      typename Block::iterator lb(block.lower_bound(j));
      typename Block::iterator it(lb);
      // no entries have ever been computed for this column, polynomial and level
      if (lb == block.end() ||
          block.key_comp()(j, lb->first))
          {
             //cout << "Punkt 3" << endl; 
//              cout << "insert a new level: " << j << ", in the block: " << lambda << endl;
              typedef typename Block::value_type value_type;
              it = block.insert(lb, value_type(j, Subblock()));


              Subblock& subblock(it->second);

              IntersectingList nus;


              intersecting_quarklets(frame(), lambda,
                                    std::max(j, frame().j0()),
                                    j == (frame().j0()-1),
                                    nus, p);

      
  
	  // do the rest of the job
	  const double d1 = problem->D(lambda);
//          cout << d1 << endl;
          
	    for (typename IntersectingList::const_iterator it2(nus.begin()), itend2(nus.end());
		 it2 != itend2; ++it2) {
                
	      const double entry = problem->a(*it2, lambda);
 //             cout << *it2 << ", " << entry << endl;
	      typedef typename Subblock::value_type value_type_subblock;
	      if (entry != 0.) {
 //                 cout << "ja3" << endl;
		subblock.insert(subblock.end(), value_type_subblock(number(*it2,jmax), entry));
	      //w.add_coefficient(*it2, (entry / (d1 * problem->D(*it2))) * factor);
 //               cout << *it2 << ": " << (*it2).number() << endl;
 //             cout <<  number(*it2,jmax) << endl;
 //               cout << w.size() << endl;
                w[number(*it2,jmax)] += (entry / (d1*problem->D(*it2))) * factor;
                //cout << "ja5" << endl;
	      }
 //             cout << "ja1" << endl;
	    }
//	  cout << "ja2" << endl;
         }

        
	
      else {
//          cout << "Einträge schon berechnet" << endl;
	// level already exists --> extract level from cache

	
          //cout << "Iterator: " << it->first << endl;
          Subblock& subblock(it->second);
        
	
	const double d1 = problem->D(lambda);
	
	//cout << "Länge w: " << w.size() << endl;
	  for (typename Subblock::const_iterator it2(subblock.begin()), itend2(subblock.end());
	       it2 != itend2; ++it2) {
//              w.add_coefficient(*(frame().get_wavelet(it2->first)),
// 			      (it2->second / (d1 * problem->D( *(frame().get_wavelet(it2->first)) )))  * factor);
                    //cout << it2->second << endl;
                    //cout << it2->first << endl;
              
                    w[it2->first] += (it2->second / (d1*problem->D(*(frame().get_quarklet(it2->first))))) * factor;
                }
	
      }// end else

      
  }


  template <class PROBLEM>
  double
  CachedQuarkletProblem<PROBLEM>::norm_A() const
  {
    if (normA == 0.0) {
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "CachedProblem()::norm_A() called..." << endl;
#endif

      std::set<Index> Lambda;
      const int j0 = frame().j0();
      const int jmax = std::min(frame().get_jmax_(),j0+2);
      const int pmax = std::min(frame().get_pmax_(),2);
      
      int p = 0;
      
      for (Index lambda = frame().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == frame().last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == frame().last_wavelet(jmax,p)){
            ++p;
            lambda = frame().first_generator(j0,p);
        }
        else
            ++lambda;
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
      
      //cout << "Eigenwerte: " << evals << endl;
      //cout << "Eigenvektoren: " << endl << evecs << endl;
      normA = evals(evals.size()-1);
      //cout << "normA: " << normA << endl;
      int i = 0;
      while(abs(evals(i))<1e-1){
          ++i;
      }
      normAinv = 1./evals(i);
      //cout << "inversnorm: " << normAinv << endl;
#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "... done!" << endl;
#endif
    }

    return normA;
  }
   
  template <class PROBLEM>
  double
  CachedQuarkletProblem<PROBLEM>::norm_Ainv() const
  {
    if (normAinv == 0.0) {
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "CachedProblem()::norm_Ainv() called..." << endl;
#endif

      std::set<Index> Lambda;
      const int j0 = frame().j0();
      const int jmax = std::min(frame().get_jmax_(),j0+2);
      const int pmax = std::min(frame().get_pmax_(),2);
//      cout << "jmax: " << jmax << endl;
//      cout << "pmax: " << pmax << endl;
      
      int p = 0;
      
      for (Index lambda = frame().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == frame().last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == frame().last_wavelet(jmax,p)){
            ++p;
            lambda = frame().first_generator(j0,p);
        }
        else
            ++lambda;
      }
      
      
//      cout << "Schritt 1" << endl;
      SparseMatrix<double> A_Lambda;
      Matrix<double> evecs;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
//      cout << "Matrix aufgestellt" << endl;
      //cout << A_Lambda << endl;
     Vector<double> evals;
      
#if 0
      double help;
      unsigned int iterations;
      cout << "Beginne Normberechnung" << endl;
      LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      
      normAinv = 1./help;
      cout << "Norm " << normA << endl;
      cout << "Inv Norm " << normAinv << endl;
      cout << iterations << endl;
#else
      //LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
      SymmEigenvalues(A_Lambda, evals, evecs);
      
      cout << "Eigenwerte: " << evals << endl;
      //cout << "Eigenvektoren: " << endl << evecs << endl;
      normA = evals(evals.size()-1);
      cout << "normA: " << normA << endl;
      int i = 0;
      while(abs(evals(i))<1e-1){
          ++i;
      }
      normAinv = 1./evals(i);
      cout << evals(i) << endl;
      cout << "inversnorm: " << normAinv << endl;
#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
      cout << "... done!" << endl;
#endif
    }

    return normAinv;
  }
  
  
}
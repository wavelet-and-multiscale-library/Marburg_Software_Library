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
  double
  CachedQuarkletProblem<PROBLEM>::a(const Index& lambda,
			    const Index& nu) const
  {
      double r = 0;
      
      WaveletBasis mybasis(basis());
      
      
      const int jmax = mybasis.get_jmax_();
      const int pmax = mybasis.get_pmax_();
      const int waveletsonplevel = ((mybasis.last_wavelet(jmax)).number()+1);
      
    
      const int lambda_num = lambda.number() + lambda.p() * waveletsonplevel;
      const int nu_num = nu.number()+ nu.p() * waveletsonplevel;
      
      //cout << "Punkt 1" << endl;
      // BE CAREFUL: KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
      typedef typename Index::type_type generator_type;
      int j = (lambda.e() == generator_type()) ? (lambda.j()-1) : lambda.j();
      
      // check wether entry has already been computed
      typedef std::list<Index> IntersectingList;

      // search for column 'nu'
      typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
      typename ColumnCache::iterator col_it(col_lb);
      
      if (col_lb == entries_cache.end() ||
	  entries_cache.key_comp()(nu_num, col_lb->first))
	
	{
	  // insert a new column
	  typedef typename ColumnCache::value_type value_type;
	  col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
	}
	 //cout << "Punkt 2" << endl; 
      Column& col(col_it->second);
      
      // check wether the level 'lambda' belongs to has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);
      // no entries have ever been computed for this column and this level
      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
         //cout << "Punkt 3" << endl; 
	  // compute whole level block
	  // #### ONLY CDD COMPRESSION STRATEGY IMPLEMENTED ####
	  // #### MAYBE WE ADD TRUNK FOR STEVENSON APPROACH ####
	  // #### LATER.                                    ####
	  
	  // insert a new level
	  typedef typename Column::value_type value_type;
	  it = col.insert(lb, value_type(j, Block()));

	  Block& block(it->second);	  

	  IntersectingList nus;
          
          for(int p=0; p<=pmax; p++){
	   
            intersecting_wavelets(basis(), nu,
                                  std::max(j, basis().j0()),
                                  j == (basis().j0()-1),
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
//                cout << (*it).number()+(*it).p()*waveletsonplevel << endl;
                typedef typename Block::value_type value_type_block;
                if (entry != 0.) {
                    block.insert(block.end(), value_type_block((*it).number()+(*it).p()*waveletsonplevel, entry));
                    if ((*it).number()+(*it).p()*waveletsonplevel == lambda_num) {
                        r = entry;
                    }
                }
            }
          }
          //cout << "nu" << nu << endl;
          //cout << block.end() << endl;
	}
      // level already exists --> extract row corresponding to 'lambda'
      else {
	  //cout << "ja" << endl;
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
        
    
    return r;
  }
  
  template <class PROBLEM>
  void
  CachedQuarkletProblem<PROBLEM>::add_level (const Index& lambda,
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
//      cout << "bin in addlevel drin" << endl;
//      cout << lambda << endl;
//      cout << j << endl;
      
      WaveletBasis mybasis(basis());
      
      const int minplevel = std::max(0, lambda.p() + 1  - (int) pow(2,(J-b*abs(j-lambda.j())/a)));
      const int maxplevel = std::min(lambda.p() - 1  + (int) pow(2,(J-b*abs(j-lambda.j())/a)), pmax);
      const int waveletsonplevel = ((mybasis.last_wavelet(jmax)).number()+1);
      const int lambda_num = lambda.number() + lambda.p() * waveletsonplevel;
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

      // check wether the level has already been calculated
      typename Column::iterator lb(col.lower_bound(j));
      typename Column::iterator it(lb);

      if (lb == col.end() ||
	  col.key_comp()(j, lb->first))
	{
	  // no entries have ever been computed for this column and this level
//         cout << "keine Einträge berechnet" << endl;
          
	  
	  // insert a new level
	  typedef typename Column::value_type value_type;
	  it = col.insert(lb, value_type(j, Block()));

	  
          Block& block(it->second);
          
          
          IntersectingList nus;
          
          
          //DKOR strategy
          
//          cout << lambda << endl;
//          cout << j << endl;
//          cout << J << endl;
          //cout << minplevel << endl;
          //cout << maxplevel << endl;
          
         for(int p = minplevel; p<=maxplevel; ++p){
             //cout << "Hallo" << endl;
          intersecting_wavelets(basis(), lambda,
				std::max(j, basis().j0()),
				j == (basis().j0()-1),
				nus, p);
          
          
          
          

	  
	  // do the rest of the job
	  const double d1 = problem->D(lambda);
          //cout << d1 << endl;
          
	    for (typename IntersectingList::const_iterator it2(nus.begin()), itend2(nus.end());
		 it2 != itend2; ++it2) {
                
	      const double entry = problem->a(*it2, lambda);
              //cout << *it2 << ", " << entry << endl;
	      typedef typename Block::value_type value_type_block;
	      if (entry != 0.) {
                  //cout << "ja3" << endl;
		block.insert(block.end(), value_type_block((*it2).number()+(*it2).p()* waveletsonplevel, entry));
	      //w.add_coefficient(*it2, (entry / (d1 * problem->D(*it2))) * factor);
                //cout << *it2 << ": " << (*it2).number() << endl;
		//cout << (*it2).number()+(*it2).p()*((1<<jmax+1)-1) << endl;
                
                w[(*it2).number()+(*it2).p()* waveletsonplevel] += (entry / (d1*problem->D(*it2))) * factor;
                //cout << "ja5" << endl;
	      }
              //cout << "ja1" << endl;
	    }
	  //cout << "ja2" << endl;
         }

        
	}
      else {
          //cout << "Einträge schon berechnet" << endl;
	// level already exists --> extract level from cache

	
          //cout << "Iterator: " << it->first << endl;
          Block& block(it->second);
        
	
	const double d1 = problem->D(lambda);
	
	//cout << "Länge w: " << w.size() << endl;
	  for (typename Block::const_iterator it2(block.begin()), itend2(block.end());
	       it2 != itend2; ++it2) {
//              w.add_coefficient(*(problem->basis().get_wavelet(it2->first)),
// 			      (it2->second / (d1 * problem->D( *(problem->basis().get_wavelet(it2->first)) )))  * factor);
                    //cout << it2->second << endl;
                    //cout << it2->first << endl;
              if(it2->first >= w.size()){
                  cout << "Feeehler" << endl;
                  cout << "level: " << j << endl;
                  cout << it->first << endl;
                  cout << col_it->first << endl;
                  //cout << *col_it << endl;
                  for (typename Block::const_iterator it3(block.begin()), itend2(block.end());
	       it3 != itend2; ++it3) {
                      
                  cout << it3->first << endl;
                  
                  }
              }
              assert(it2->first<w.size());
                    w[it2->first] += (it2->second / (d1*problem->D(*(problem->basis().get_wavelet(it2->first))))) * factor;
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
      const int j0 = basis().j0();
      const int jmax = j0+2;
      const int pmax = std::min(basis().get_pmax_(),2);
      //const int pmax = 0;
      int p = 0;
      
      for (Index lambda = problem->basis().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == problem->basis().last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == problem->basis().last_wavelet(jmax,p)){
            ++p;
            lambda = problem->basis().first_generator(j0,p);
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
      const int j0 = basis().j0();
      const int jmax = j0+2;
      int pmax = std::min(basis().get_pmax_(),2);
      
      int p = 0;
      
      for (Index lambda = problem->basis().first_generator(j0,0);;) {
	Lambda.insert(lambda);
	if (lambda == problem->basis().last_wavelet(jmax,pmax)) break;
        //if (i==7) break;
        if (lambda == problem->basis().last_wavelet(jmax,p)){
            ++p;
            lambda = problem->basis().first_generator(j0,p);
        }
        else
            ++lambda;
      }
      
      
      //cout << "Schritt 1" << endl;
      SparseMatrix<double> A_Lambda;
      Matrix<double> evecs;
      setup_stiffness_matrix(*this, Lambda, A_Lambda);
      //cout << "Matrix aufgestellt" << endl;
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

    return normAinv;
  }
}
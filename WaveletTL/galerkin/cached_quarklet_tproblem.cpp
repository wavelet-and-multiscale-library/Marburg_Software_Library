// implementation for cached_quarklet_tproblem.h

namespace WaveletTL
{
    template <class PROBLEM>
    CachedQuarkletTProblem<PROBLEM>::CachedQuarkletTProblem(PROBLEM* P,
                                            const double estnormA,
                                            const double estnormAinv)
    : problem(P), normA(estnormA), normAinv(estnormAinv)
    {
    }
    
    template <class PROBLEM>
    double
    CachedQuarkletTProblem<PROBLEM>::a(const Index& la,
			       const Index& mu) const
    {
//        printf("enter a");
        double r = 0;

        if (problem->local_operator())
        {
            const Index* lambda= &la;
            const Index* nu = &mu;
            
            const int lambda_num = lambda->number();
//            cout << endl << "Lambda: " << lambda  << endl;
            const int nu_num = nu->number();
//            cout << "Nu: " << nu << ", " << endl;
//            const unsigned int pmax = problem->frame().get_pmax();          

            // Be careful, there is no generator level in the cache! The situation from the MRA setting:
            // KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
            // does not hold in the tensor setting. Generators and quarklets on
            // the minimal level are thrown together in one index set (componentwise),
            // this also applies to tensors of generators on the lowest level

//            level_type lambda_key/*,first_level(frame().j0())*/;
//            level_type lambda_key;
//            polynomial_type lambda_p(lambda.p());
//            polynomial_type lambda_p5(lambda.p());
//            polynomial_type lambda_p1(lambda.p());
//            polynomial_type lambda_p2(lambda.p());
//            polynomial_type lambda_p3(lambda.p());
//            polynomial_type lambda_p4(lambda.p());
//            for (int k=0;k<2;k++)
//            {
//                lambda_key[k] = lambda->j()[k]-3/*frame().j0()[k]*//*first_level[k]*/;
//            }

            int blocknumber((lambda->p()).number());
//            int blocknumber(0);
//HACK: due to performance reasons we compute the relativelevelnumber directly instead of calling the number() routine
            int jmulti = lambda->j()[0] +lambda->j()[1]-2*frame().j0()[0];
            int relativelevelnumber = jmulti*(jmulti+1)/2 + lambda->j()[0]-frame().j0()[0];
//            if(lambda_key.number()==relativelevelnumber)
//                cout << "ja" << *lambda << endl;
//            else{
//                cout << "nein" << *lambda << endl;
//                abort();
//            }
            
            int subblocknumber(relativelevelnumber);
//            int subblocknumber(0);

            // check wether entry has already been computed
//            typedef std::list<Index> IntersectingList;
//            typedef std::list<Index*> IntersectingPointerList;
            typedef std::list<int> IntersectingList;

            // search for column 'mu'
            typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
            typename ColumnCache::iterator col_it(col_lb);

            if (col_lb == entries_cache.end() ||
                entries_cache.key_comp()(nu_num, col_lb->first))
            {
//                cout << "Neue Spalte " << nu << endl;
                // insert a new column
                typedef typename ColumnCache::value_type value_type;

#pragma omp critical //(cache_insert)
                {
                    col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
                }

            }
            

            Column& col(col_it->second);

            // check wether the polynomial 'lambda' belongs to has already been calculated
            typename Column::iterator block_lb(col.lower_bound(blocknumber));
            typename Column::iterator block_it(block_lb);

            if (block_lb == col.end() ||
                col.key_comp()(blocknumber, block_lb->first))
            {
                
                
//                cout << "Neuer block" << endl;
                typedef typename Column::value_type value_type;

#pragma omp critical //(column_insert)
                {
                    block_it = col.insert(block_lb, value_type(blocknumber, Block()));
                }

            }
            
            Block& block(block_it->second);
            
            
            
            
            
            // check wether the level 'lambda' belongs to has already been calculated
            typename Block::iterator lb(block.lower_bound(subblocknumber));
            typename Block::iterator it(lb);

            if (lb == block.end() ||
                block.key_comp()(subblocknumber, lb->first))
            {
//                cout << "subblocknumber: " << subblocknumber << ", " << lambda << endl;
            
                // no entries have ever been computed for this column, polynomial and level
                // compute whole level block
                // insert a new level
//                cout << "Neuer subblock" << endl;
//                cout << "Lambda: " << lambda  << endl;
//                cout << "Nu: " << nu << ", " << endl << endl;

                typedef typename Block::value_type value_type;

#pragma omp critical //(block_insert)
                {
                    it = block.insert(lb, value_type(subblocknumber, Subblock()));
                }

                
                Subblock& subblock(it->second);

//                
                // there are no Generators
//                IntersectingList nus;
//                IntersectingPointerList nus;
                IntersectingList nus;
//                cout << 
//                    cout << lambda->j() << endl;
                intersecting_quarklets(frame(),
                                       *nu,
                                       //std::max(j, frame().j0()),
                                       lambda->j(),
                                       nus,
                                       lambda->p());
                
//                intersecting_quarklets(frame(),
//                                       nu_num,
//                                       //std::max(j, frame().j0()),
//                                       lambda->j(),
//                                       nus,
//                                       lambda->p());
                
//                frame();
//                frame();
//                frame();
//                frame();
//                                       nu_num;
//                                       //std::max(j, frame().j0()),
//                                       lambda->j();
//                                       nus;
//                                       lambda->p();
                
                // compute entries
//                cout << "bin hier"<<endl;
                //cout<<nus.size()<<endl;

                
//                for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)
//                 for (typename IntersectingPointerList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)   
                double entry;
#if PARALLEL_A==1
#pragma omp parallel num_threads(NUM_THREADS) 
#endif
                {
#if PARALLEL_A==1
#pragma omp for private(entry) schedule(dynamic) nowait
                 for(int n=0;n<nus.size();n++){
                  typename IntersectingList::const_iterator it(nus.begin());
                  advance(it,n);              
#else
                for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)    
                {
#endif

                    //hier second compression einbauen 
                    //if(!(strategy==tensor_second) || (strategy==tensor_second && (dist<radius/2 || intersect_singular_support(frame(),*(frame().get_quarklet(*it)), *nu)) ) ) {
                    //cout<<*it<<endl;
//                    const double entry = 0;
//                    const double entry = problem->a(*it, *nu);
//                    const double entry = problem->a(*it, *nu);
//                    const double entry = problem->a(*(frame().get_quarklet((*it).number())), *nu);
                      entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                    const double entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                    *(frame().get_quarklet(*it));


                    typedef typename Subblock::value_type value_type_subblock;
                    if (fabs(entry) > 1e-16 ) 
//                            if (entry != 0.)
                    {

//                        subblock.insert(subblock.end(), value_type_subblock((*it).number(), entry));

#pragma omp critical //(subblock_insert)
                        {
                        subblock.insert(subblock.end(), value_type_subblock(*it, entry));
                        }

//                                cout << *it << endl;
//                        if ((int)(*it).number() == lambda_num)
                        if (*it == lambda_num)    
                        {
                            r = entry;
                        }
                    }
                    //} //second compression ende
                }
            }
            }
            // level already exists --> extract row corresponding to 'lambda'
            else
//#pragma omp critical
            {
//                cout << "Ergebnis aus L2 Cache"  << endl;
                
                Subblock& subblock(it->second);

                //typename Block::iterator block_lb(block.lower_bound(lambda));
                typename Subblock::iterator subblock_lb(subblock.lower_bound(lambda_num));
                typename Subblock::iterator subblock_it(subblock_lb);
                // level exists, but in row 'lambda' no entry is available ==> entry must be zero
                if (subblock_lb == subblock.end() ||
                    subblock.key_comp()(lambda_num, subblock_lb->first))
                  {
                    r = 0;
//                    cout << "Ergebnis ist 0" << endl;
                  }
                else {
                    
                  r = subblock_it->second;
//                  cout << "Ergenis: " << r << endl;
                }
            }
        }
//       cout<<"exit a"<<endl;  
      return r;
      
    }
     
    //version for second compression 
        //mu is fixed, la is variable
    template <class PROBLEM>
    double
    CachedQuarkletTProblem<PROBLEM>::a(const Index& la,
			       const Index& mu, const int radius, const CompressionStrategy strategy, const double A, const double B) const
    {
//        printf("bin hier \n");
        double r = 0;

        if (problem->local_operator())
        {
            const Index* lambda= &la;
            const Index* nu = &mu;
            
            double dist=A*log(1+abs(nu->p()[0]-lambda->p()[0])+abs(nu->p()[1]-lambda->p()[1]))+B*(abs(nu->j()[0]-lambda->j()[0])+abs(nu->j()[1]-lambda->j()[1]));
//          cout<<"dist="<<dist<<"  radius="<<radius<<endl;
            
            const int lambda_num = lambda->number();
//            cout << endl << "Lambda: " << lambda  << endl;
            const int nu_num = nu->number();
//            cout << "Nu: " << nu << ", " << endl;
//            const unsigned int pmax = problem->frame().get_pmax();          

            // Be careful, there is no generator level in the cache! The situation from the MRA setting:
            // KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
            // does not hold in the tensor setting. Generators and quarklets on
            // the minimal level are thrown together in one index set (componentwise),
            // this also applies to tensors of generators on the lowest level

//            level_type lambda_key/*,first_level(frame().j0())*/;
//            level_type lambda_key;
//            polynomial_type lambda_p(lambda.p());
//            polynomial_type lambda_p5(lambda.p());
//            polynomial_type lambda_p1(lambda.p());
//            polynomial_type lambda_p2(lambda.p());
//            polynomial_type lambda_p3(lambda.p());
//            polynomial_type lambda_p4(lambda.p());
//            for (int k=0;k<2;k++)
//            {
//                lambda_key[k] = lambda->j()[k]-3/*frame().j0()[k]*//*first_level[k]*/;
//            }

            int blocknumber((lambda->p()).number());
//            int blocknumber(0);
//HACK: due to performance reasons we compute the relativelevelnumber directly instead of calling the number() routine
            int jmulti = lambda->j()[0] +lambda->j()[1]-2*frame().j0()[0];
            int relativelevelnumber = jmulti*(jmulti+1)/2 + lambda->j()[0]-frame().j0()[0];
//            if(lambda_key.number()==relativelevelnumber)
//                cout << "ja" << *lambda << endl;
//            else{
//                cout << "nein" << *lambda << endl;
//                abort();
//            }
            
            int subblocknumber(relativelevelnumber);
//            int subblocknumber(0);

            // check wether entry has already been computed
//            typedef std::list<Index> IntersectingList;
//            typedef std::list<Index*> IntersectingPointerList;
            typedef std::list<int> IntersectingList;

            // search for column 'mu'
            typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
            typename ColumnCache::iterator col_it(col_lb);

            if (col_lb == entries_cache.end() ||
                entries_cache.key_comp()(nu_num, col_lb->first))
            {
//                cout << "Neue Spalte " << nu << endl;
                // insert a new column
                typedef typename ColumnCache::value_type value_type;
                
#pragma omp critical //(cache_insert)
                {
                    col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
                }

            }
            

            Column& col(col_it->second);

            // check wether the polynomial 'lambda' belongs to has already been calculated
            typename Column::iterator block_lb(col.lower_bound(blocknumber));
            typename Column::iterator block_it(block_lb);

            if (block_lb == col.end() ||
                col.key_comp()(blocknumber, block_lb->first))
            {
                
                
//                cout << "Neuer block" << endl;
                typedef typename Column::value_type value_type;

#pragma omp critical //(columen_insert)
                {
                    block_it = col.insert(block_lb, value_type(blocknumber, Block()));
                }

            }
            
            Block& block(block_it->second);
            
            
            
            
            
            // check wether the level 'lambda' belongs to has already been calculated
            typename Block::iterator lb(block.lower_bound(subblocknumber));
            typename Block::iterator it(lb);

            if (lb == block.end() ||
                block.key_comp()(subblocknumber, lb->first))
            {
//                cout << "subblocknumber: " << subblocknumber << ", " << lambda << endl;
            
                // no entries have ever been computed for this column, polynomial and level
                // compute whole level block
                // insert a new level
//                cout << "Neuer subblock" << endl;
//                cout << "Lambda: " << lambda  << endl;
//                cout << "Nu: " << nu << ", " << endl << endl;

                typedef typename Block::value_type value_type;

#pragma omp critical //(block_insert)
                {
                    it = block.insert(lb, value_type(subblocknumber, Subblock()));
                }

                
                Subblock& subblock(it->second);

//                
                // there are no Generators
//                IntersectingList nus;
//                IntersectingPointerList nus;
                IntersectingList nus;
//                cout << 
//                    cout << lambda->j() << endl;
                intersecting_quarklets(frame(),
                                       *nu,
                                       //std::max(j, frame().j0()),
                                       lambda->j(),
                                       nus,
                                       lambda->p());
                
//                intersecting_quarklets(frame(),
//                                       nu_num,
//                                       //std::max(j, frame().j0()),
//                                       lambda->j(),
//                                       nus,
//                                       lambda->p());
                
//                frame();
//                frame();
//                frame();
//                frame();
//                                       nu_num;
//                                       //std::max(j, frame().j0()),
//                                       lambda->j();
//                                       nus;
//                                       lambda->p();
                
                // compute entries
                //cout << "bin hier"<<endl;
                //cout<<nus.size()<<endl;

                
//                for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)
//                 for (typename IntersectingPointerList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)  
                double entry;
#if PARALLEL_A==1
#pragma omp parallel num_threads(NUM_THREADS) 
#endif
                {
#if PARALLEL_A==1
#pragma omp for private(entry) schedule(dynamic) nowait
                 for(int n=0;n<nus.size();n++){
                  typename IntersectingList::const_iterator it(nus.begin());
                  advance(it,n);              
#else
                for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)    
                {
#endif
                    //hier second compression einbauen 
                    if(!(strategy==tensor_second) || (strategy==tensor_second && (dist<=radius*0.5 || intersect_singular_support(frame(),*(frame().get_quarklet(*it)), *nu)) ) ) 
                    {
                    //cout<<*it<<endl;
//                    const double entry = 0;
//                    const double entry = problem->a(*it, *nu);
//                    const double entry = problem->a(*it, *nu);
//                    const double entry = problem->a(*(frame().get_quarklet((*it).number())), *nu);
                      entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                    const double entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                    *(frame().get_quarklet(*it));


                    typedef typename Subblock::value_type value_type_subblock;
                    if (fabs(entry) > 1e-16 ) 
//                            if (entry != 0.)
                    {

//                        subblock.insert(subblock.end(), value_type_subblock((*it).number(), entry));

#pragma omp critical //(subblock_insert)
                        {
                        subblock.insert(subblock.end(), value_type_subblock(*it, entry));
                        }

//                                cout << *it << ", " << (*it).number() <<endl;
//                        if ((int)(*it).number() == lambda_num)
                        if (*it == lambda_num)    
                        {
                            r = entry;
                        }
                    }
                    } //second compression ende
//                    else if(strategy==tensor_second && dist>radius/2 && !intersect_singular_support(frame(),*(frame().get_quarklet(*it)), *nu)){
//                        cout<<"second compression: dropping additional entries"<<endl;
//                        cout<<*nu<<"  ,  "<<*(frame().get_quarklet(*it))<<endl;
//                    }
                }
                }
            }
            // level already exists --> extract row corresponding to 'lambda'
            else
            {
//                cout << "Ergebnis aus L2 Cache"  << endl;
             
                Subblock& subblock(it->second);

                //typename Block::iterator block_lb(block.lower_bound(lambda));
                typename Subblock::iterator subblock_lb(subblock.lower_bound(lambda_num));
                typename Subblock::iterator subblock_it(subblock_lb);
                // level exists, but in row 'lambda' no entry is available ==> entry must be zero
                if (subblock_lb == subblock.end() ||
                    subblock.key_comp()(lambda_num, subblock_lb->first))
                  {
                    r = 0;
//                    cout << "Ergebnis ist 0" << endl;
                  }
                else {

                  r = subblock_it->second;
                    
//                  cout << "Ergenis: " << r << endl;
                }
                
            }
        }

      return r;
      
    }    

//    template <class PROBLEM>
//    double
//    CachedQuarkletTProblem<PROBLEM>::a(const Index& la,
//			       const Index& mu) const
//    {
//        
//        double r = 0;
//
//        if (problem->local_operator())
//        {
//            const Index* lambda= &la;
//            const Index* nu = &mu;
//            
//            const int lambda_num = lambda->number();
//            const int nu_num = nu->number();
////            const unsigned int pmax = problem->frame().get_pmax();          
//
//            // Be careful, there is no generator level in the cache! The situation from the MRA setting:
//            // KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
//            // does not hold in the tensor setting. Generators and quarklets on
//            // the minimal level are thrown together in one index set (componentwise),
//            // this also applies to tensors of generators on the lowest level
//
//            level_type lambda_key(lambda->j()),first_level(frame().j0());
//            polynomial_type lambda_p(lambda->p());
//            for (int k=0;k<space_dimension;k++)
//            {
//                lambda_key[k] = lambda_key[k]-first_level[k];
//            }
////TODO (PERFORMANCE): store the numbers of all levels up to jmax, do not compute anything here:
//            int blocknumber((lambda->p()).number());
//            int subblocknumber(lambda_key.number());
//
//            // check wether entry has already been computed
//            typedef std::list<int> IntersectingList;
//
//            // search for column 'mu'
//            typename ColumnCache::iterator col_lb(entries_cache.lower_bound(nu_num));
//            typename ColumnCache::iterator col_it(col_lb);
//
//            if (col_lb == entries_cache.end() ||
//                entries_cache.key_comp()(nu_num, col_lb->first))
//            {
////                cout << "Neue Spalte " << nu << endl;
//                // insert a new column
//                typedef typename ColumnCache::value_type value_type;
//                col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
//            }
//            
//
//            Column& col(col_it->second);
//
//            // check wether the polynomial 'lambda' belongs to has already been calculated
//            typename Column::iterator block_lb(col.lower_bound(blocknumber));
//            typename Column::iterator block_it(block_lb);
//
//            if (block_lb == col.end() ||
//                col.key_comp()(blocknumber, block_lb->first))
//            {
//                
//                
//
//                typedef typename Column::value_type value_type;
//                block_it = col.insert(block_lb, value_type(blocknumber, Block()));
//            }
//            
//            Block& block(block_it->second);
//            
//            
//            
//            
//            
//            // check wether the level 'lambda' belongs to has already been calculated
//            typename Block::iterator lb(block.lower_bound(subblocknumber));
//            typename Block::iterator it(lb);
//
//            if (lb == block.end() ||
//                block.key_comp()(subblocknumber, lb->first))
//            {
////                cout << "subblocknumber: " << subblocknumber << ", " << lambda << endl;
//            
//                // no entries have ever been computed for this column, polynomial and level
//                // compute whole level block
//                // insert a new level
//                
//
//                typedef typename Block::value_type value_type;
//                it = block.insert(lb, value_type(subblocknumber, Subblock()));
//                Subblock& subblock(it->second);
//
//                if ((subblocknumber == 0) && (lambda->j() == first_level))
//                {
//                    // Insert Generators and Quarklets
//
//                    // also add Generators to the Cache
//                    // this case has to be considered seperatly because
//                    // intersecting_quarklets divides the case of a frame
//                    // function made entirely of generators and the cache does not.
//
//                    // TODO : produce nicer looking code:
//                    IntersectingList nusG,nusW;
//                        intersecting_quarklets(frame(), *nu,
//                                               lambda->j(),
////                                          true, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tframe_support!
//                                          nusG,
//                                          lambda->p());
//                        intersecting_quarklets(frame(), *nu,
//                                          lambda->j(),
////                                          false, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tframe_support!
//                                          nusW,
//                                          lambda->p());
//                        // compute entries
//                        for (typename IntersectingList::const_iterator it(nusG.begin()), itend(nusG.end());it != itend; ++it)
//                        {
////                            const double entry = problem->a(*it, nu);
//                            const double entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                            typedef typename Subblock::value_type value_type_subblock;
//                            if (fabs(entry) > 1e-16 ) 
////                            if (entry != 0.)
//                            {
//                                subblock.insert(subblock.end(), value_type_subblock(*it, entry));
////                                if ((int)(*it).number() == lambda_num)
////                                {
////                                    r = entry;
////                                }
//                                if (*it == lambda_num)    
//                                {
//                                    r = entry;
//                                }
//                            }
//                        }
//                        for (typename IntersectingList::const_iterator it(nusW.begin()), itend(nusW.end());it != itend; ++it)
//                        {
////                            const double entry = problem->a(*it, nu);
//                            const double entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                            typedef typename Subblock::value_type value_type_subblock;
//                            if (fabs(entry) > 1e-16 ) 
////                            if (entry != 0.)
//                            {
//                                // Insertion should be efficient, since the quarklets coma after the generators
//                                subblock.insert(subblock.end(), value_type_subblock(*it, entry));
////                                if ((int)(*it).number() == lambda_num)
////                                {
////                                    r = entry;
////                                }
//                                if (*it == lambda_num)    
//                                {
//                                    r = entry;
//                                }
//                            }
//                        }
//                    }
//                
//                    else
//                    {
//                        // there are no Generators
//                        IntersectingList nus;
//                        intersecting_quarklets(frame(),
//                                               *nu,
//                                               //std::max(j, frame().j0()),
//                                               lambda->j(),
////                                               false, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tframe_support!
//                         //                    == (frame().j0()-1),
//                                               nus,
//                                               lambda->p());
//                    // compute entries
//                        for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)
//                        {
////                            const double entry = problem->a(*it, nu);
//                            const double entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                            
//                            
//                            typedef typename Subblock::value_type value_type_subblock;
//                            if (fabs(entry) > 1e-16 ) 
//    //                            if (entry != 0.)
//                            {
//                                subblock.insert(subblock.end(), value_type_subblock(*it, entry));
////                                if ((int)(*it).number() == lambda_num)
////                                {
////                                    r = entry;
////                                }
//                                if (*it == lambda_num)    
//                                {
//                                    r = entry;
//                                }
//                            }
//                        }
//                    }
//            }
//            // level already exists --> extract row corresponding to 'lambda'
//            else
//            {
//                
//                
//                Subblock& subblock(it->second);
//
//                //typename Block::iterator block_lb(block.lower_bound(lambda));
//                typename Subblock::iterator subblock_lb(subblock.lower_bound(lambda_num));
//                typename Subblock::iterator subblock_it(subblock_lb);
//                // level exists, but in row 'lambda' no entry is available ==> entry must be zero
//                if (subblock_lb == subblock.end() ||
//                    subblock.key_comp()(lambda_num, subblock_lb->first))
//                  {
//                    r = 0;
//                  }
//                else {
//                  r = subblock_it->second;
//                }
//            }
//        }
//         // end of local operator
//
//        
//      return r;
//    }

    template <class PROBLEM>
    void
    CachedQuarkletTProblem<PROBLEM>::add_ball(const Index& lambda,
				      //InfiniteVector<double, Index>& w,
				      Vector<double>& w,
				      const int radius,
				      const double factor,
				      const int maxlevel,
                                      const CompressionStrategy strategy,
                                      const bool precond,
                                      const int maxpolynomial,
                                      const double A,
                                      const double B) const
    {
        /*
         * For dim=1,2 it is easy to describe all levels in the (1-norm) ball of radius around lambda.
         * For higher dimensions this is done recursivly.
         */
        
//        cout<<"add_ball wird ausgefuehrt"<<endl;
//        cout<<lambda<<endl;
        
        assert (space_dimension <= 2);// we only support dimensions up to 2
        
        double d1 = precond?D(lambda):1.0;
        polynomial_type p;
#if 0        
        if (space_dimension == 1)
        {
            
        
            int j0 = this->frame().j0()[0];
            
            Index mu;
            const int minlevel = max (j0, lambda.j()[0]-(int) (radius / B));
            const int finalmaxlevel = min(lambda.j()[0]+ (int) (radius / B), maxlevel);
            for (int level = minlevel; level  < finalmaxlevel ;level++)
            {
                int minpolynomial = std::max(0, lambda.p()[0] + 1  - (int) pow(2,(radius-B*abs(level-lambda.j()[0])/A)));
                int finalmaxpolynomial = std::min(lambda.p()[0] - 1  + (int) pow(2,(radius-B*abs(level-lambda.j()[0])/A)), maxpolynomial);
                
                
                p[0]=minpolynomial;
                for(;p[0]<=finalmaxpolynomial; ++p){
                
                    // @ hier weitermachen: for schleife für p's einfügen
                    mu = this->frame().first_quarklet(level, p);
                    // the result of the following call is that the cache holds the whole level block corresponding to all quarklets with level = |mu|
                    // ideally one should replace this method call with a slightly modified version of the code from a(,)
                    // using the first generator on level j0 would mean no difference, since block 0 contains quarklets and generators
                    a(mu,lambda);
                    
                    // add the level
                    Subblock& subblock (entries_cache[lambda.number()][p[0]][level-j0]);
    #if 1
                    for (typename Subblock::const_iterator it(subblock.begin()), itend(subblock.end()); it != itend; ++it)
                    {



                        // high caching strategy:
                        // The call of D fills the diagonal block of each column in block
                        w[it->first]=w[it->first]+factor*(it->second)/( precond? (d1*D(this->frame().get_quarklet(it->first))):1.0);
                        // low caching strategy:
                        // by calling the uncached Version of D no column but the one of lambda is altered. This Version is slower than the high caching strategy.
                        //w[it->first]=w[it->first]+factor*(it->second)/d1/(problem->D(this->frame().get_quarklet(it->first)));
                    }
                }
#else
                // mid caching strategy:
                // maybe it is too much effort to compute the whole diagonal block (as in the high caching strategy).
                // For the columns in block other than the one of lambda we use D_2 to only compute and cache a(mu,mu) (a single entry)
                // Do do this a bool variable has to be introduced marking that only one entry in mus column is stored.
                // this caching strategy is useful for sparse signals with active coefficients on high levels
                // ... needs to be tested though ...
                typename Block::const_iterator it(block.begin()), itend(block.end());
                w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->frame().get_quarklet(it->first));
                ++it;
                for (; it != itend; ++it)
                {
                    w[it->first]=w[it->first]+factor*(it->second)/d1/D_already(this->frame().get_quarklet(it->first));
                }
#endif
            }
        }
#endif
        if (space_dimension == 2)
        {

            // The ball can be described of levellines consisting of levels with the same multidegree
            // The first is determined with the distance of lambda.j and j0. The last with maxlevel
            // The first level in a levelline is determined with minx = min(j0[0], lambda.j[0]-radius)
            // The last level in a levelline is determined with miny = min(j0[1], lambda.j[1]-radius)

            
//            Index mu;
            level_type j0(this->frame().j0());
            
            int lambdaline = multi_degree(lambda.j());
            int actualradius =(int)(radius/B);
//            int lambdaline = lambda.j()[0]+lambda.j()[1];
            int lowestline = multi_degree(j0);
//            int lowestline = j0[0]+j0[1];
            int dist2j0=lambdaline-lowestline;
            int dist2maxlevel=maxlevel-lambdaline;
            level_type currentlevel, leveldiffj0;
            polynomial_type currentpolynomial;
            
            int xstart,xend,ystart,subblocknumber;
            // iterate the levellines. offset relative to lambdas levelline

            {                     
            for (int offset = -std::min(dist2j0,actualradius); offset <= std::min(dist2maxlevel,actualradius); offset++)
            {
                // iterate over the levels on the levelline
                xstart = lambda.j()[0]-actualradius+ceil((actualradius+offset)/2.0); //x coordinate of the first level, ignoring restrictions by j0[0]
                xend = lambda.j()[0]+floor((actualradius+offset)/2.0); // same for the last
                ystart = lambda.j()[1]+floor((actualradius+offset)/2.0); // and for the second dimension
                // the first level in the ball is denoted with steps=0.
                // The first level on the current levelline in the ball may have steps >0.
                for (int steps = max(0,j0[0]-xstart); steps <= min(xend-xstart,ystart-j0[1]) ; steps++ )
                {
                    currentlevel[0]= xstart+steps;
                    leveldiffj0[0]=currentlevel[0]-frame().j0()[0];
                    currentlevel[1]= ystart-steps;
                    leveldiffj0[1]=currentlevel[1]-frame().j0()[1];
                    subblocknumber = currentlevel[0]-j0[0] + ((lambdaline-lowestline+offset)*(lambdaline-lowestline+offset+1))/2;
                    
                    
//                    int dist2p0=multi_degree(lambda.p());
//                    int dist2maxpolynomial=maxlevel-dist2p0;
                    int jxdiff = currentlevel[0]-lambda.j()[0];
                    int jydiff = currentlevel[1]-lambda.j()[1];
                    int gap = abs(jxdiff) + abs(jydiff); 
                    double boundary;// = pow(2,(radius-B*gap/A))/((1+lambda.p()[0])*(1+lambda.p()[1])); //ÄNDERN 1 bezieht sich auf (1+p_1')*(1+p_2')
                    if(A==B){
                        boundary = 1<<(radius-gap)/((1+lambda.p()[0])*(1+lambda.p()[1]));
                    }
                    else{
                        boundary = pow(2,(radius-B*gap/A))/((1+lambda.p()[0])*(1+lambda.p()[1])); //ÄNDERN 1 bezieht sich auf (1+p_1')*(1+p_2')
                    } 
                    //                    cout << "boundary: " << boundary << endl;
//                    cout << "maxpolynomial: " << maxpolynomial << endl;
//                    cout << "gap: " << gap << endl;
//                    cout << "radius: " << radius << endl;
                    
                    for (currentpolynomial[0]=0; currentpolynomial[0]<= std::min((int)boundary-1,maxpolynomial);currentpolynomial[0]++){
                        
                        for (currentpolynomial[1]=0; currentpolynomial[1]<= std::min((int)(boundary/(1+currentpolynomial[0]))-1,maxpolynomial-currentpolynomial[0]);currentpolynomial[1]++){
 //                           Index ind=frame().get_quarklet(frame().get_first_wavelet_numbers()[leveldiffj0.number()]+currentpolynomial.number() * frame().get_Nablasize());
//                            double dist=A*log(1+abs(ind.p()[0]-lambda.p()[0])+abs(ind.p()[1]-lambda.p()[1]))+B*(abs(ind.j()[0]-lambda.j()[0])+abs(ind.j()[1]-lambda.j()[1]));
//                            cout<<"dist="<<dist<<endl;
//                            double dist2=A*log(1+abs(currentpolynomial[0]-lambda.p()[0])+abs(currentpolynomial[1]-lambda.p()[1]))+B*(abs(currentlevel[0]-lambda.j()[0])+abs(currentlevel[1]-lambda.j()[1]));
//                            cout<<"dist2="<<dist2<<endl;
                            //here second compression is applied
                            //if(!(strategy==tensor_second) || (strategy==tensor_second && (dist2<radius/2 || intersect_singular_support(frame(),lambda,*(frame().get_quarklet(frame().get_first_wavelet_numbers()[leveldiffj0.number()]+currentpolynomial.number() * frame().get_Nablasize())) ) )) ){
//                            mu = this->frame().first_quarklet(currentlevel, currentpolynomial);
//                            a(mu,lambda);
//                            cout<<"computing a"<<endl;
                            if(strategy==tensor_second){
                               a(*(frame().get_quarklet(frame().get_first_wavelet_numbers()[leveldiffj0.number()]+currentpolynomial.number() * frame().get_Nablasize())),lambda, radius, strategy, A, B); 
                            }
                            else{
                            a(*(frame().get_quarklet(frame().get_first_wavelet_numbers()[leveldiffj0.number()]+currentpolynomial.number() * frame().get_Nablasize())),lambda);
                            }
//                            cout<<"nach a"<<endl;
                            //                            cout << mu << endl;
//                            int dist2p0 = multi_degree(lambda.p());
//                            int blocknumber = currentpolynomial[0] + ((dist2p0+offset)*(dist2p0+offset+1))/2;;
                            


                            // add the level
//#pragma omp critical
//                            {
                            Subblock& subblock (entries_cache[lambda.number()][currentpolynomial.number()][subblocknumber]);
                            
//                            cout<<"nach subblock"<<endl;
        #if 1
                            for (typename Subblock::const_iterator it(subblock.begin()), itend(subblock.end()); it != itend; ++it)
                            {
                                // high caching strategy
                                const double d2=D(*(frame().get_quarklet(it->first)));
                                
#pragma omp critical //(write_output)
                                {
                                w[it->first]=w[it->first]+factor*(it->second)/( precond? (d1*d2):1.0);
                                }
//                                w[it->first]=w[it->first]+factor*(it->second)/( precond? (d1*D(this->frame().get_quarklet(it->first))):1.0);

                                // low caching strategy
                                //w[it->first]=w[it->first]+factor*(it->second)/d1/(problem->D(this->frame().get_quarklet(it->first)));
                            }
//                            }
        #else
                            
                            // mid caching strategy
                            typename Block::const_iterator it(block.begin()), itend(block.end());
                            w[it->first]=w[it->first]+factor*(it->second)/d1/D(this->frame().get_quarklet(it->first));
                            ++it;
                            for (; it != itend; ++it)
                            {
                                w[it->first]=w[it->first]+factor*(it->second)/d1/D_already(this->frame().get_quarklet(it->first));
                            }
#endif
                    //}//end of second compression
                    }
                    }
                }
            }

        }      

        }
//        else // dim > 2. iteration over all levels in 'range' is done recursivly for arbitrary dimensions
//        {
//            add_level_recurse(lambda,w,radius,factor,lambda.j(),0,maxlevel,true,d1,precond);
//        }
    }



    template <class PROBLEM>
    double
    CachedQuarkletTProblem<PROBLEM>::norm_A() const
    {
        if (normA == 0.0)
        {
            //const int offset = 2;
            int offsetj, offsetp;
            switch (space_dimension)
            {
                case 1:
                    offsetj = 2, offsetp = std::min((int)problem->frame().get_pmax(),2);
                    break;
                case 2:
                    offsetj = 1, offsetp = std::min((int)problem->frame().get_pmax(),2);
                    break;
                default:
                    offsetj = 0, offsetp = 0;
            }
        
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "CachedQuarkletTProblem()::norm_A() called..." << endl;
#endif

            set<Index> Lambda;
            polynomial_type p;
            //cout << "cached_tproblem.norm_A :: last quarklet = " << (problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offset)) << endl;
            for (Index lambda ( problem->frame().first_generator(problem->frame().j0(), p) );;)
            {
                
//                cout << "index: " << lambda << endl;
                Lambda.insert(lambda);
                if (lambda == problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offsetj, p) ){
                    ++p;
                    if ((int)multi_degree(p)>offsetp) break;
                    lambda = problem->frame().first_generator(problem->frame().j0(), p);
                }
                else
                    ++lambda;
                   
            }
//            for (typename set<Index>::const_iterator it(Lambda.begin()), itend(Lambda.end());it != itend; ++it)
//                        {   
//                
//                cout << "Lambda Index: " << *it << endl;
//                            
//                        }
            
            
            SparseMatrix<double> A_Lambda;
            cout << "begin setup_stiffness_matrix" << endl;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
            cout << "end setup_stiffness_matrix" << endl;
            
//#if 1
            
            
            A_Lambda.compress(1e-10);
//            unsigned int iterations;
//            double help;
//            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            Matrix<double> evecs;
            Vector<double> evals;
            SymmEigenvalues(A_Lambda, evals, evecs);
            int i = 0;
            while(abs(evals(i))<1e-2){
                ++i;
            }
            normA = evals(evals.size()-1);
            normAinv = 1./evals(i);
            
//#else
//          Vector<double> xk(Lambda.size(), false);
//          xk = 1;
//          unsigned int iterations;
//          normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
//#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "... done!" << endl;
#endif
        }
        return normA;
    }

    template <class PROBLEM>
    void
    CachedQuarkletTProblem<PROBLEM>::normtest(unsigned int offsetj, unsigned int offsetp) const
    {
        //if (normA == 0.0)
        {
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "CachedQuarkletTProblem()::norm_A() called..." << endl;
#endif
            set<Index> Lambda;
            //cout << "cached_tproblem.norm_A :: last quarklet = " << (problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offset)) << endl;
            polynomial_type p;
            //cout << "cached_tproblem.norm_A :: last quarklet = " << (problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offset)) << endl;
            for (Index lambda ( problem->frame().first_generator(problem->frame().j0(), p) );;)
            {
                Lambda.insert(lambda);
                if (lambda == problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offsetj, p) ){
                    ++p;
                    if (multi_degree(p)>offsetp) break;
                    lambda = problem->frame().first_generator(problem->frame().j0(), p);
                }
                else
                    ++lambda;
                   
            }
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
//#if 1
//            double help;
//            unsigned int iterations;
//            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            A_Lambda.compress(1e-10);
//            unsigned int iterations;
//            double help;
//            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            Matrix<double> evecs;
            Vector<double> evals;
            SymmEigenvalues(A_Lambda, evals, evecs);
            int i = 0;
            while(abs(evals(i))<1e-2){
                ++i;
            }
            normA = evals(evals.size()-1);
            normAinv = 1./evals(i);
            
//#else
//          Vector<double> xk(Lambda.size(), false);
//          xk = 1;
//          unsigned int iterations;
//          normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
//#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "... done!" << endl;
#endif
        }
        cout << "normtest:: offsetj = " << offsetj <<  " offsetp = " << offsetp << "  normA = " << normA << " normAinv = " << normAinv << " Kondition = " << (normA*normAinv) << endl;
    }

    template <class PROBLEM>
    double
    CachedQuarkletTProblem<PROBLEM>::norm_Ainv() const
    {
        //cout << " this is norm_Ainv(). Verbosity is " << _WAVELETTL_CACHEDPROBLEM_VERBOSITY << endl;
        if (normAinv == 0.0)
        {
            int offsetj, offsetp;
            switch (space_dimension)
            {
                case 1:
                    offsetj = 2, offsetp = std::min((int)problem->frame().get_pmax(),2);
                    break;
                case 2:
                    offsetj = 1, offsetp = std::min((int)problem->frame().get_pmax(),2);
                    break;
                default:
                    offsetj = 0, offsetp = 0;
            }
            
            
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "CachedQuarkletTProblem()::norm_Ainv() called..." << endl;
            cout << "verbosity is "<< _WAVELETTL_CACHEDPROBLEM_VERBOSITY << endl;
#endif
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >=2
            Index tempindex(problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offset));
            cout << "number of last quarklet is "<<tempindex.number()<<endl;
#endif
            set<Index> Lambda;
            polynomial_type p;
            //cout << "cached_tproblem.norm_A :: last quarklet = " << (problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offset)) << endl;
            for (Index lambda ( problem->frame().first_generator(problem->frame().j0(), p) );;)
            {
                Lambda.insert(lambda);                
                if (lambda == problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offsetj, p) ){
                    ++p;
                    if ((int)multi_degree(p)>offsetp) break;
                    lambda = problem->frame().first_generator(problem->frame().j0(), p);
                }
                else
                    ++lambda;
            }
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
//#if 1
//            double help;
//            unsigned int iterations;
//            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
//            normAinv = 1./help;
            
            A_Lambda.compress(1e-10);
//            unsigned int iterations;
//            double help;
//            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            Matrix<double> evecs;
            Vector<double> evals;
            SymmEigenvalues(A_Lambda, evals, evecs);
            int i = 0;
            while(abs(evals(i))<1e-2){
                ++i;
            }
            normA = evals(evals.size()-1);
            normAinv = 1./evals(i);
//#else
//            Vector<double> xk(Lambda.size(), false);
//            xk = 1;
//            unsigned int iterations;
//            normAinv = InversePowerIteration(A_Lambda, xk, 1e-6, 200, iterations);
//#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "... done!" << endl;
#endif
        }
        return normAinv;
    }

    template <class PROBLEM>
    void
    CachedQuarkletTProblem<PROBLEM>::set_f(const Function<PROBLEM::space_dimension>* fnew)
    {
        problem->set_f(fnew);
    }



    
}

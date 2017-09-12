// implementation for cached_quarklet_ldomain_problem.h

namespace WaveletTL
{
    template <class PROBLEM>
    CachedQuarkletLDomainProblem<PROBLEM>::CachedQuarkletLDomainProblem(PROBLEM* P,
                                            const double estnormA,
                                            const double estnormAinv)
    : problem(P), normA(estnormA), normAinv(estnormAinv)
    {
    }

    template <class PROBLEM>
    double
    CachedQuarkletLDomainProblem<PROBLEM>::a(const Index& la,
			       const Index& mu) const
    {
        
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
                col_it = entries_cache.insert(col_lb, value_type(nu_num, Column()));
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
                block_it = col.insert(block_lb, value_type(blocknumber, Block()));
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
                it = block.insert(lb, value_type(subblocknumber, Subblock()));
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
#if 0
#pragma omp parallel for shared(subblock)
                for(int i=0;i<nus.size();i++){
                    typename IntersectingList::const_iterator it(nus.begin());
                    advance(it,i);
#else
                
//                for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)
//                 for (typename IntersectingPointerList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)   
                for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)    
                {
#endif
                    //cout<<*it<<endl;
//                    const double entry = 0;
//                    const double entry = problem->a(*it, *nu);
//                    const double entry = problem->a(*it, *nu);
//                    const double entry = problem->a(*(frame().get_quarklet((*it).number())), *nu);
                        const double entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                    const double entry = problem->a(*(frame().get_quarklet(*it)), *nu);
//                    *(frame().get_quarklet(*it));


                    typedef typename Subblock::value_type value_type_subblock;
                    if (fabs(entry) > 1e-16 ) 
//                            if (entry != 0.)
                    {
//                        subblock.insert(subblock.end(), value_type_subblock((*it).number(), entry));
                        subblock.insert(subblock.end(), value_type_subblock(*it, entry));
//                                cout << *it << ", " << (*it).number() <<endl;
//                        if ((int)(*it).number() == lambda_num)
                        if (*it == lambda_num)    
                        {
                            r = entry;
                        }
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

    template <class PROBLEM>
    void
    CachedQuarkletLDomainProblem<PROBLEM>::add_ball(const Index& lambda,
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
        //cout <<"entering add_ball"<<endl; 
        //cout << lambda << " : " << radius << " : " << factor << " : " << maxlevel<<endl;
        
        
        assert (space_dimension <= 2);// we only support dimensions up to 2
        
        double d1 = precond?D(lambda):1.0;
        polynomial_type p;
        
        
        

        // The ball can be described of levellines consisting of levels with the same multidegree
        // The first is determined with the distance of lambda.j and j0. The last with maxlevel
        // The first level in a levelline is determined with minx = min(j0[0], lambda.j[0]-radius)
        // The last level in a levelline is determined with miny = min(j0[1], lambda.j[1]-radius)


        //Index mu;
        level_type j0(this->frame().j0());

        int lambdaline = multi_degree(lambda.j());
        int actualradius =(int)(radius/B);
//            int lambdaline = lambda.j()[0]+lambda.j()[1];
        int lowestline = multi_degree(j0);
//            int lowestline = j0[0]+j0[1];
        int dist2j0=lambdaline-lowestline;
        int dist2maxlevel=maxlevel-lambdaline;
        level_type currentlevel;
        polynomial_type currentpolynomial;

        int xstart,xend,ystart,subblocknumber;
        // iterate the levellines. offset relative to lambdas levelline
#if 0
#pragma omp parallel for
#endif
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
                currentlevel[1]= ystart-steps;
                subblocknumber = currentlevel[0]-j0[0] + ((lambdaline-lowestline+offset)*(lambdaline-lowestline+offset+1))/2;


//                    int dist2p0=multi_degree(lambda.p());
//                    int dist2maxpolynomial=maxlevel-dist2p0;
                int jxdiff = currentlevel[0]-lambda.j()[0];
                int jydiff = currentlevel[1]-lambda.j()[1];
                int gap = abs(jxdiff) + abs(jydiff);
                double boundary;
                if(A==B){
                    boundary = 1<<(radius-gap)/((1+lambda.p()[0])*(1+lambda.p()[1]));
                }
                else{
                    boundary = pow(2,(radius-B*gap/A))/((1+lambda.p()[0])*(1+lambda.p()[1])); //ÄNDERN 1 bezieht sich auf (1+p_1')*(1+p_2')
                }
                    //cout << "boundary: " << boundary << endl;
//                    cout << "maxpolynomial: " << maxpolynomial << endl;
//                    cout << "gap: " << gap << endl;
//                    cout << "radius: " << radius << endl;

                for (currentpolynomial[0]=0; currentpolynomial[0]<= std::min((int)boundary-1,maxpolynomial);currentpolynomial[0]++)
                {

                    for (currentpolynomial[1]=0; currentpolynomial[1]<= std::min((int)(boundary/(1+currentpolynomial[0]))-1,maxpolynomial-currentpolynomial[0]);currentpolynomial[1]++)
                    {

                        //mu = this->frame().first_quarklet(currentlevel, currentpolynomial);
                        // the result of the following call is that the cache holds the whole level subblock corresponding to all 
                        // quarklets with level and polynomial degree of mu
                        a(this->frame().first_quarklet(currentlevel, currentpolynomial),lambda);
//                            cout << mu << endl;
//                            int dist2p0 = multi_degree(lambda.p());
//                            int blocknumber = currentpolynomial[0] + ((dist2p0+offset)*(dist2p0+offset+1))/2;;



                        // add the level
                        Subblock& subblock (entries_cache[lambda.number()][currentpolynomial.number()][subblocknumber]);
    
                        for (typename Subblock::const_iterator it(subblock.begin()), itend(subblock.end()); it != itend; ++it)
                        {
                            // high caching strategy
                            w[it->first]=w[it->first]+factor*(it->second)/( precond? (d1*D(this->frame().get_quarklet(it->first))):1.0);

                            // low caching strategy
                            //w[it->first]=w[it->first]+factor*(it->second)/d1/(problem->D(this->frame().get_quarklet(it->first)));
                        }
    
                    }
                }
            }
        }
        
    }
    



    

    


    template <class PROBLEM>
    double
    CachedQuarkletLDomainProblem<PROBLEM>::norm_A() const
    {
        if (normA == 0.0)
        {
            //const int offset = 2;
            int offsetj, offsetp;
            
            offsetj = std::min(1,problem->frame().get_jmax()-(int)multi_degree(problem->frame().j0())); 
            offsetp = std::min(problem->frame().get_pmax(),2);
                    
        
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "CachedQuarkletLDomainProblem()::norm_A() called..." << endl;
#endif

            set<Index> Lambda;
            polynomial_type p;
//            cout << "Hier" << endl;
//            cout << "offsetj: " << offsetj << endl;
//            cout << problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offsetj,p) << endl;
//            cout << "cached_ldomain_problem.norm_A :: last quarklet = " << (problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offsetj,p)) 
//                    << (problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offsetj,p)).number() << endl;
            for (Index lambda ( problem->frame().first_generator(problem->frame().j0(), p) );;)
            {
                
//                cout << "index: " << lambda << endl;
                Lambda.insert(lambda);
//                cout << "hier" <<endl;
                problem->frame().last_quarklet(multi_degree(problem->frame().j0())+offsetj, p);
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
//                cout << "Lambda Index: " << *it << ", " << (*it).number() << endl;
//                            
//                        }
            
            
            SparseMatrix<double> A_Lambda;
            cout << "begin setup_stiffness_matrix" << endl;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
            cout << "end setup_stiffness_matrix, computing eigenvalues" << endl;
            
#if 0
            
            
            A_Lambda.compress(1e-10);
//            unsigned int iterations;
//            double help;
 //           LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            Matrix<double> evecs;
            Vector<double> evals;
            SymmEigenvalues(A_Lambda, evals, evecs);
            int i = 0;
            while(abs(evals(i))<1e-2){
                ++i;
            }
            normA = evals(evals.size()-1);
            normAinv = 1./evals(i);
            
#else
          Vector<double> xk(Lambda.size(), false);
          xk = 1;
          unsigned int iterations;
//          normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
         double lambdamax=PowerIteration<Vector<double>,SparseMatrix<double> >(A_Lambda, xk, 1e-3, 100, iterations);
         double lambdamin=InversePowerIteration<Vector<double>,SparseMatrix<double> >(A_Lambda, xk, 1e-2,  1e-3, 100, iterations);
         normA=lambdamax;
         normAinv=1./lambdamin;
#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "... done!" << endl;
#endif
        }
        return normA;
    }

    template <class PROBLEM>
    void
    CachedQuarkletLDomainProblem<PROBLEM>::normtest(unsigned int offsetj, unsigned int offsetp) const
    {
        //if (normA == 0.0)
        
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "CachedQuarkletLDomainProblem()::norm_A() called..." << endl;
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
            double lambdamin, lambdamax;
//            if (offsetp>0){
//    //            double lambdamin;
//    //            unsigned int iterations;
//    //            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
//                A_Lambda.compress(1e-10);
//    //            unsigned int iterations;
//    //            double help;
//    //            LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
//                Matrix<double> evecs;
//                Vector<double> evals;
//                SymmEigenvalues(A_Lambda, evals, evecs);
//    //            cout << "Eigenwerte: " << evals << endl;
//                int i = 0;
//                while(abs(evals(i))<1e-2){
//                    ++i;
//                }
//                normA = evals(evals.size()-1);
//                lambdamin = evals(i);
//                normAinv = 1./lambdamin;
////#else
////        unsigned int iterations;
////        LanczosIteration(A_Lambda, 1e-6, lambdamin, lambdamax, 200, iterations);
////        normAinv = 1./lambdamin;
////        normA= lambdamax;
////#endif
//            
//            }
            
                Vector<double> xk(Lambda.size(), false);
                Vector<double> yk(Lambda.size(), false);
                xk = 1, yk = 1;
                unsigned int iterations;
        //          normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
                lambdamax=PowerIteration<Vector<double>,SparseMatrix<double> >(A_Lambda, xk, 1e-3, 100, iterations);

                lambdamin=InversePowerIteration<Vector<double>,SparseMatrix<double> >(A_Lambda, yk, 1e-2,  1e-3, 100, iterations);
                normA=lambdamax;
                normAinv=1./lambdamin;
            
           
//#else
//          Vector<double> xk(Lambda.size(), false);
//          xk = 1;
//          unsigned int iterations;
//          normA = PowerIteration(A_Lambda, xk, 1e-6, 100, iterations);
//#endif

#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "... done!" << endl;
#endif
        
        cout << "normtest:: offsetj = " << offsetj <<  " offsetp = " << offsetp << "  normA = " << normA << " normAinv = " << normAinv << endl
                << " kleinster Eigenwert = " << lambdamin <<" Kondition = " << (normA*normAinv) << endl;
    }

    template <class PROBLEM>
    double
    CachedQuarkletLDomainProblem<PROBLEM>::norm_Ainv() const
    {
        //cout << " this is norm_Ainv(). Verbosity is " << _WAVELETTL_CACHEDPROBLEM_VERBOSITY << endl;
        if (normAinv == 0.0)
        {
            int offsetj, offsetp;
            offsetj = std::min(1,problem->frame().get_jmax()-(int)multi_degree(problem->frame().j0()));;
            offsetp = std::min((int)problem->frame().get_pmax(),2);
                    
            
            
#if _WAVELETTL_CACHEDPROBLEM_VERBOSITY >= 1
            cout << "CachedQuarkletLDomainProblem()::norm_Ainv() called..." << endl;
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
    CachedQuarkletLDomainProblem<PROBLEM>::set_f(const Function<PROBLEM::space_dimension>* fnew)
    {
        problem->set_f(fnew);
    }
    
}

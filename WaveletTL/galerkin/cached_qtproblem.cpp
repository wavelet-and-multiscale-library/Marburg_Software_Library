// implementation for cached_qtproblem.h



namespace WaveletTL
{
      
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::CachedQTProblem(QTBASIS* basis,
                //const Array1D<FixedMatrix<double, DER_ONEDIMHAARCOUNT> > & agencoeffs,
                //const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > & qgencoeffs,
                const double normA,
                const double normAinv)
    : basis_(basis), normA_(normA), normAinv_(normAinv)
    {
        list<int> temp_list;
        for (unsigned int p=0; p< basis_->get_nop(); p++)
        {
            temp_list.push_back(p);
        }
        
        for (unsigned int i=0; i<DIM; ++i)
        {
            temp_list.sort(index_cmp< Array1D<MultiIndex<int,DIM> > > (basis_->j0(),i));
            int p_= temp_list.front();
            int temp_int = basis_->j0()[p_][i];
            patch_with_minimal_j0_[i].clear();
            for (list<int>::const_iterator it(temp_list.begin()), itend(temp_list.end()); it != itend; ++it)
            {
                if (temp_int == basis_->j0()[*it][i] )
                {
                    patch_with_minimal_j0_[i].push_back(*it);
                }
                else
                {
                    break;
                }
            }
            assert (patch_with_minimal_j0_[i].size() > 0);
        }
        if ( (normA_ == 0.0) || (normAinv_ == 0.0))
        {
            std::set<int> Lambda;
            const int offset(min (basis_->get_jmax() - multi_degree(basis_->j0()[0]), (unsigned int)0)); // offset of 2 to kep computation effort low)
            
            MultiIndex<int, DIM> temp_jmax(basis_->j0()[0]);
            temp_jmax[0]=temp_jmax[0]+offset; 
            
            unsigned int num(basis_->get_levelnum(temp_jmax, basis_->get_nop()-1));
            Index temp_ind(basis_->last_wavelet(num));
            num = temp_ind.number()+1;
            for (unsigned int n=0;n<num; ++n)
            {
                Lambda.insert(n);
            }
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
            double help;
            unsigned int iterations;
            //LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            LanczosIteration(A_Lambda, 1e-6, help, normA_, 200, iterations);
            normAinv_ = 1./help;
        }
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::CachedQTProblem(QTBASIS* basis,
                const Array1D<FixedMatrix<double, DER_ONEDIMHAARCOUNT> > & agencoeffs,
                const Array1D<FixedMatrix<double, ONEDIMHAARCOUNT> > & qgencoeffs,
                const Function<DIM>* f,
                const char* rhs_filename,
                const double normA,
                const double normAinv)
    : basis_(basis), agencoeffs_(agencoeffs), qgencoeffs_(qgencoeffs), f_(f), normA_(normA), normAinv_(normAinv)
    {
        compute_D();
        if (f != 0)
        {
            compute_rhs(rhs_filename);
        }
        else
        {
            load_rhs(rhs_filename);
        }
        list<int> temp_list;
        for (unsigned int p=0; p< basis_->get_nop(); p++)
        {
            temp_list.push_back(p);
        }
        
        for (unsigned int i=0; i<DIM; ++i)
        {
            temp_list.sort(index_cmp< Array1D<MultiIndex<int,DIM> > > (basis_->j0(),i));
            int p_= temp_list.front();
            int temp_int = basis_->j0()[p_][i];
            patch_with_minimal_j0_[i].clear();
            for (list<int>::const_iterator it(temp_list.begin()), itend(temp_list.end()); it != itend; ++it)
            {
                if (temp_int == basis_->j0()[*it][i] )
                {
                    patch_with_minimal_j0_[i].push_back(*it);
                }
                else
                {
                    break;
                }
            }
            assert (patch_with_minimal_j0_[i].size() > 0);
        }
        if ( (normA_ == 0.0) || (normAinv_ == 0.0))
        {
            std::set<int> Lambda;
            const int offset(min (basis_->get_jmax() - multi_degree(basis_->j0()[0]), (unsigned int)0)); // offset of 2 to kep computation effort low)
            
            MultiIndex<int, DIM> temp_jmax(basis_->j0()[0]);
            temp_jmax[0]=temp_jmax[0]+offset; 
            
            unsigned int num(basis_->get_levelnum(temp_jmax, basis_->get_nop()-1));
            Index temp_ind(basis_->last_wavelet(num));
            num = temp_ind.number()+1;
            for (unsigned int n=0;n<num; ++n)
            {
                Lambda.insert(n);
            }
            SparseMatrix<double> A_Lambda;
            setup_stiffness_matrix(*this, Lambda, A_Lambda);
            double help;
            unsigned int iterations;
            //LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
            LanczosIteration(A_Lambda, 1e-6, help, normA_, 200, iterations);
            normAinv_ = 1./help;
        }
    };

    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::compute_D()
    {
        diagonal_cache_.resize(basis_->degrees_of_freedom());
        for (unsigned int i = 0; i< basis_->degrees_of_freedom();i++)
        {
            diagonal_cache_[i] = sqrt(a(i,i));
        }
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::compute_rhs(const char* rhs_filename)
    {
        // precompute the right-hand side on a fine level
        InfiniteVector<double,int> fhelp;
        fcoeffs_precond_unsorted_.resize(0);
        fcoeffs_precond_unsorted_.resize(basis_->degrees_of_freedom());
        f_precond_norm_sqr_ = 0;
        double coeff;
        for (unsigned int i = 0; i< basis_->degrees_of_freedom();i++)
        {
            coeff = basis_->integrate(f_,i)/D(i);
            if (fabs(coeff)>1e-15)
            {
                fhelp.set_coefficient(i, coeff);
                fcoeffs_precond_unsorted_[i] = coeff;
                f_precond_norm_sqr_ += coeff*coeff;
            }
        }
        // sort the coefficients into fcoeffs
        fcoeffs_precond_sorted_.resize(0); // clear eventual old values
        fcoeffs_precond_sorted_.resize(fhelp.size());
        unsigned int id(0);
        for (typename InfiniteVector<double,int>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                it != itend; ++it, ++id)
        {
            fcoeffs_precond_sorted_[id] = std::pair<int,double>(it.index(), *it);
        }
        sort(fcoeffs_precond_sorted_.begin(), fcoeffs_precond_sorted_.end(), typename InfiniteVector<double,int>::decreasing_order());
        // store on disc
        std::ofstream ofs(rhs_filename,std::ofstream::binary);
        writeVectorToFile(fcoeffs_precond_unsorted_, ofs);
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::load_rhs(const char* rhs_filename) //)
    {
        std::ifstream ifs(rhs_filename, std::ifstream::binary);
        readVectorFromFile(fcoeffs_precond_unsorted_,ifs);
        assert (basis_->degrees_of_freedom() == fcoeffs_precond_unsorted_.size());
        InfiniteVector<double,int> fhelp;
        f_precond_norm_sqr_ = 0;
        double coeff;
        for (unsigned int i = 0; i< basis_->degrees_of_freedom();i++)
        {
            //coeff = basis_->integrate(f_,i)/D(i);
            coeff = fcoeffs_precond_unsorted_[i];
            if (coeff != 0)
            {
                fhelp.set_coefficient(i, coeff);
                //fcoeffs_precond_unsorted_[i] = coeff;
                f_precond_norm_sqr_ += coeff*coeff;
            }
        }
        
        // sort the coefficients into fcoeffs
        fcoeffs_precond_sorted_.resize(0); // clear eventual old values
        fcoeffs_precond_sorted_.resize(fhelp.size());
        unsigned int id(0);
        for (typename InfiniteVector<double,int>::const_iterator it(fhelp.begin()), itend(fhelp.end());
                it != itend; ++it, ++id)
        {
            fcoeffs_precond_sorted_[id] = std::pair<int,double>(it.index(), *it);
        }
        sort(fcoeffs_precond_sorted_.begin(), fcoeffs_precond_sorted_.end(), typename InfiniteVector<double,int>::decreasing_order());
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    double 
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::a(const int munum,
            const int nunum)
    {
        
        //für festes nu, berechne alle mu2 mit selbem Level wie mu die nu überlappen
        double r = 0;
        Index mu(basis_->get_wavelet(munum));
        Index nu (basis_->get_wavelet(nunum));

        const MultiIndex<int, DIM> nu_j(nu.j()), nu_e(nu.e()), nu_k(nu.k()),mu_j(mu.j()), mu_e(mu.e()), mu_k(mu.k());
        const unsigned int nu_p(nu.p()), mu_p(mu.p());

        MultiIndex<unsigned int, DIM> intinfo;
        //bool temp_b(basis_->get_LMR_info(nunum,mu_j,mu_p, intinfo));
        bool temp_b(basis_->get_LMR_info(nunum,mu_p, intinfo));
        if (!temp_b) 
            return 0;
        FixedArray1D<entries,DIM> integralshares; // store the value of the one dimensional integrals mu_i nu_i against Haar generators
        int kmingeni, kmaxgeni, kminwavi, kmaxwavi;
        unsigned int nui_basisnum, mui_basisnum;
        bool mu_min_type_i;
        for (unsigned int i=0; i<DIM; i++)
        {
            nui_basisnum = (((basis_->get_bc()[nu_p][2*i])?0:2) + ((basis_->get_bc()[nu_p][2*i+1])?0:1));
            mui_basisnum = (((basis_->get_bc()[mu_p][2*i])?0:2) + ((basis_->get_bc()[mu_p][2*i+1])?0:1));
            mu_min_type_i = (mu_j[i]==basis_->j0()[mu_p][i]) ? true:false;
            ColumnCache* cachePointer;
            Column* cachePointerGen;
            switch (intinfo[i])
            {
                case 0:
                case 4:
                case 8: // LL MM RR -> type 1
                    cachePointer = & typeIcache_[4*nui_basisnum + mui_basisnum];
                    if (mu_min_type_i)
                        cachePointerGen = & typeIcachedGens_[4*nui_basisnum + mui_basisnum];
                    break;
                case 1:
                case 7: 
                case 2:
                case 6: // LM RM LR RL -> type 2
                    cachePointer = &typeIIcache_[4*(nui_basisnum-1) + mui_basisnum];
                    if (mu_min_type_i)
                        cachePointerGen = & typeIIcachedGens_[4*(nui_basisnum-1) + mui_basisnum];
                    break;
                case 3: // ML -> type 3             
                    // EXTREME CAUTION:: mui_basisnum - (intinfo[i] == 5)?0:1 produces garbage!
                    assert (mui_basisnum != 1);
                    cachePointer = &typeIIIcache_[4*nui_basisnum + mui_basisnum - 1];
                    if (mu_min_type_i)
                        cachePointerGen = & typeIIIcachedGens_[4*nui_basisnum + mui_basisnum - 1];
                    break;
                case 5: // MR -> type 3                
                    assert (mui_basisnum != 2);
                    cachePointer = &typeIIIcache_[4*nui_basisnum + ((mui_basisnum == 1)?0:3)];
                    if (mu_min_type_i)
                        cachePointerGen = & typeIIIcachedGens_[4*nui_basisnum + ((mui_basisnum == 1)?0:3)];
                    break;
                default:
                    abort();
            }
            // search for column 'mu_i'
            typename QTBASIS::IntervalBasis::Index nui(nu_j[i],nu_e[i],nu_k[i],basis_->get_bases_infact()[nui_basisnum]);
            unsigned int nui_num = nui.number();
            typename ColumnCache::iterator col_lb(cachePointer->lower_bound(nui_num));
            typename ColumnCache::iterator col_it(col_lb);
            if (col_lb == cachePointer->end() ||
                cachePointer->key_comp()(nui_num, col_lb->first))
            {
                // the nui-th column has not been requested so far
                // insert a new column and continue with the blocks
                typedef typename ColumnCache::value_type value_type;
                col_it = cachePointer->insert(col_lb, value_type(nui_num, Column()));
            }
            // col_it points to the column of \nu_i
            Column& col(col_it->second);
            // check whether the level 'mu_i' belongs to has already been calculated
            int mui_levelnum (mu_j[i] - basis_->j0()[mu_p][i]);
            typename Column::iterator lb(col.lower_bound(mui_levelnum));
            typename Column::iterator it(lb);
            if (lb == col.end() ||
                col.key_comp()(mui_levelnum, lb->first))
            {
                // no entries have ever been computed for this column and this level
                // insert a new block, 
                // and then compute the whole level block 
                //      (\int psi_mui psi_nui)_mui, (\int psi_mui' psi_nui')_mui
                // for all mui that intersect nui. 
                // if mui is on the minimal level we also need to add a new Block() to
                // the generator integral cache, i.e., we need to 
                // integrate against all generators on the lowest level
                typedef typename Column::value_type value_type;
                it = col.insert(lb, value_type(mui_levelnum, Block()));
                Block& block(it->second);
                bool gen_intersection_i, wav_intersection_i;
                basis_->get_onedim_intersections(intinfo[i],
                        nu_j[i],
                        nu_e[i],
                        nu_k[i],
                        nui_basisnum,
                        (mui_levelnum == 0),
                        mu_j[i],
                        mui_basisnum,
                        kmingeni,
                        kmaxgeni,
                        kminwavi,
                        kmaxwavi,
                        gen_intersection_i,
                        wav_intersection_i);
// cleanup                
                if ((nu_j[i] == 3) && (nu_e[i] == 0 && ((nu_k[i] == 3) && ((mu_j[i] == 3) && ((nui_basisnum == 0) && (mui_basisnum == 3) && (intinfo[i] == 3))))))
                {
                    cout << "// CACHED_QTPROBLEM:: //" << endl;
                    cout << "intinfoi = " << intinfo[i] << endl;
                    cout << "kmingeni = " << kmingeni << "; "
                        << "kmaxgeni = " << kmaxgeni << "; "
                        << "kminwavi = " << kminwavi << "; "
                        << "kmaxwavi = " << kmaxwavi << "; "
                        << "gen_intersection_i = " << gen_intersection_i << "; "
                        << "wav_intersection_i = " << wav_intersection_i << endl;
                }
                if ((nu_j[i] == 3) && (nu_e[i] == 0 && ((nu_k[i] == 3) && ((mu_j[i] == 3) && ((nui_basisnum == 0) && (mui_basisnum == 3) && (intinfo[i] == 5))))))
                {
                    cout << "// CACHED_QTPROBLEM:: //" << endl;
                    cout << "intinfoi = " << intinfo[i] << endl;
                    cout << "kmingeni = " << kmingeni << "; "
                        << "kmaxgeni = " << kmaxgeni << "; "
                        << "kminwavi = " << kminwavi << "; "
                        << "kmaxwavi = " << kmaxwavi << "; "
                        << "gen_intersection_i = " << gen_intersection_i << "; "
                        << "wav_intersection_i = " << wav_intersection_i << endl;
                }
                FixedArray1D<double,ONEDIMHAARCOUNT> gram;
                FixedArray1D<double,DER_ONEDIMHAARCOUNT> der;
                if (mu_min_type_i)
                {
                    typename Column::iterator lb2(cachePointerGen->lower_bound(nui_num));
                    typename Column::iterator it2(lb2);
                    if (lb2 == cachePointerGen->end() ||
                        cachePointerGen->key_comp()(nui_num, lb2->first))
                    {
                        // no entries have ever been computed for this column and this level
                        // insert a new block, 
                        // and then compute the whole level block 
                        //      (\int psi_nui psi_mui)_nui, (\int psi_nui' psi_mui')_nui
                        // for all nui that intersect mui. 
                        // if mui is on the minimal level we also need to add a new Block() to
                        // the generator integral cache, i.e., we need to 
                        // integrate against all generators on the lowest level
                        typedef typename Column::value_type value_type;
                        it2 = cachePointerGen->insert(lb2, value_type(nui_num, Block()));
                        Block& block2(it2->second);
                        if (gen_intersection_i)
                        {
                            
                            for (int kgen = kmingeni; kgen <= kmaxgeni; ++kgen)
                            {
                                if ( (kgen == mu_k[i]) && (mu_e[i] == 0))
                                {
                                    compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                            mu_j[i], 0, kgen, mui_basisnum,
                                            gram,
                                            der);
                                    integralshares[i] = make_pair(gram,der);
                                    typedef typename Block::value_type value_type_block;
                                    block2.insert(block2.end(), value_type_block(kgen, integralshares[i]));
                                    temp_b = false;
//                                    cout << "inserted new generator Blocks. entry for mu is nontrivial" << endl;
                                }
                                else
                                {
                                    compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                            mu_j[i], 0, kgen, mui_basisnum,
                                            gram,
                                            der);
                                    
                                    typedef typename Block::value_type value_type_block;
                                    block2.insert(block2.end(), value_type_block(kgen, make_pair(gram,der)));
                                }
                            }
// CLEANUP                            
//                            if (temp_b)
//                            {
//                                cout << "inserted new generator Blocks. entry for mu is trivial. method should return 0" << endl;
//                            }
                        }
                    }
                    else
                    {
                        cout << "just inserted a wavelet block, but gen block was already existent?!!" << endl;
                        abort();
                    }
                }
                // code invariant: there exists a Block() (maybe empty) corresponding to nu and mu
                
                // if (no intersections) => current block (wavCache) remains empty
                // if we are on the lowest possible level for mui, we add an empty Block to the genCache
                // advantage: genCache[nuinum] exists and contains a block. 
                // code invariant: to access an entry: load the block and check if there is an entry (simple! same whether there is an entry or not!)
                // it would be possible to simply leave the genCache as it is.
                // advantage: nothing to do at this point
                // access to an entry: check whether there is a column in genCache. If not: value 0. If there is: get Block and take value from it.
                // This leads to 2 checks per call to a() instead of one. 
                // This is cheaper than method 1 above if nu does not intersect at all with basis functions mu.
                // However, most nus will have some intersection? In this case this variant leads to 1 additional check
                // so ... hopefully this makes the average access time to the cache faster
// cleanup                
                assert (wav_intersection_i || (!gen_intersection_i));
                if (!wav_intersection_i) return 0;
                
                // wav_intersection_i == true guarantees kminwavi <=kmaxwavi and that the values are meaningful
                for (int kwav = kminwavi; kwav <= kmaxwavi; ++kwav)
                {
                    if ( (kwav == mu_k[i]) && (mu_e[i] == 1))
                    {
                        // nu reflected?  = !(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) )
                        compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), 
                                nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                mu_j[i], 1, kwav, mui_basisnum,
                                gram,
                                der);
                        integralshares[i] = make_pair (gram,der);
                        typedef typename Block::value_type value_type_block;
                        block.insert(block.end(), value_type_block(kwav, integralshares[i]));
                        temp_b = false;
//                        cout << "inserted new wavelet Blocks. entry for mu is nontrivial" << endl;
                    }
                    else
                    {
                        compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), nu_j[i], nu_e[i], nu_k[i], nui_basisnum,
                                mu_j[i], 1, kwav, mui_basisnum,
                                gram,
                                der);
                        
                        typedef typename Block::value_type value_type_block;
                        block.insert(block.end(), value_type_block(kwav, make_pair (gram,der) ));
                    }
                }
                if (temp_b) 
                {
                    return 0;
                }
                temp_b = true;
            }
            else // column and levelblock already exist
            {
                if (mu_e[i] == 0)
                {
// CLEANUP
                    assert (mu_min_type_i);
                    typename Column::iterator lb2(cachePointerGen->lower_bound(nui_num));
                    Block& block(lb2->second);
                    typename Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    typename Block::iterator block_it(block_lb);
                    // intersections with generators exist, but 'mui_k' entry is not available ==> entry must be zero
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        return 0;
                    }
                    else 
                    {
                        integralshares[i] = block_it->second;
                    }
                }
                else
                {
                    Block& block(it->second);
                    typename Block::iterator block_lb(block.lower_bound(mu_k[i]));
                    typename Block::iterator block_it(block_lb);
                    // level exists, but in row 'mui_k' no entry is available ==> entry must be zero
                    if (block_lb == block.end() ||
                            block.key_comp()(mu_k[i], block_lb->first))
                    {
                        return 0;
                    }
                    else 
                    {
                        integralshares[i] = block_it->second;
                    }
                }
            }
        } // end of loop over dim
        // sum over all patches
        // include information about mu into LMR info:
        for (unsigned int i=0; i<DIM; ++i)
        {
            // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
            switch (intinfo[i])
            {
                case 0:
                    // LL: check whether mu is a left boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == basis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                          ||
                          ( (mu_e[i] == 1) && (mu_k[i] < basis_->get_numofbw() ) ) ) )
                    {
                        // nu is not extended left
                        intinfo[i] = 4; // MM
                    }
                    break;
                case 2:
                    // LR: check whether mu is a right boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == basis_-> get_bases(mu_p,i)->DeltaRmax(  basis_->get_j0(mu_p,i) )) ) 
                      ||
                      ( (mu_e[i] == 1) && (mu_k[i] > (1<<mu_j[i])- 1 - basis_->get_numofbw()) ) ))
                    {
                        intinfo[i] = 1; // LM
                    }
                    break;
                case 3:
                    // ML : check whether mu is a left boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == basis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                          ||
                          ( (mu_e[i] == 1) && (mu_k[i] < basis_->get_numofbw()) ) ) )
                    {
                        // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                        abort();
                    }
                    break;
                case 5:
                    // MR: check whether mu is a right boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) )) ) 
                      ||
                      ( (mu_e[i] == 1) && (mu_k[i] > (1<<mu_j[i])- 1 - basis_->get_numofbw()) ) ))
                    {
                        // nu is not extended right. There is no intersection of mu and nu! This should not happen at this point
                        cout << "cached_qtproblem::a ; i = " << i << "; intinfo = " << intinfo << endl;
                        cout << "nu_j = " << nu_j << "; nu_e = " << nu_e << "; nu_k = " << nu_k << "; nu_p = " << nu_p << endl;
                        cout << "mu_j = " << mu_j << "; mu_e = " << mu_e << "; mu_k = " << mu_k << "; mu_p = " << mu_p << endl;
                        cout << "basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) ) = " << basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) ) << endl;
                        cout << "(1<<mu_j[i])- 1 - basis_->get_numofbw() = " << (1<<mu_j[i])- 1 - basis_->get_numofbw() << endl;
                        MultiIndex<unsigned int, DIM> intinfo2;
                        
                        //temp_b = (basis_->get_LMR_info(nunum,mu_j,mu_p, intinfo2));
                        temp_b = (basis_->get_LMR_info(nunum,mu_p, intinfo2));
                        cout << "nunum = " << nunum << "; mu_j = " << mu_j << "; mu_p = " << mu_p << "; intinfo2 = " << intinfo2 << "; temp_b = " << temp_b << endl;
                        temp_b = (basis_->get_LMR_info(nunum,mu_p, intinfo2));
                        //a(munum,nunum);
                        abort();
                    }
                    break;
                case 6:
                    // RL : check whether mu is a left boundary gen/wav
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == basis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                          ||
                          ( (mu_e[i] == 1) && (mu_k[i] < basis_->get_numofbw()) ) ) )
                    {
                        // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                        intinfo[i] = 7; // RM
                    }
                    break;
                case 8:
                    // RR: check whether mu is a right boundary gen/wav
// CLEANUP
                    /*
                    if (munum == 151)
                        if (nunum == 80)
                        {
                            cout << "basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) ) = " << basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) ) << endl;
                            cout << "basis_->get_numofbw() = " << basis_->get_numofbw() << endl;
                            cout << "(1<<mu_j[i])- 1 - basis_->get_numofbw() = " << (1<<mu_j[i])- 1 - basis_->get_numofbw() << endl;
                        }
                     */
                    if (!( ( (mu_e[i] == 0) && (mu_k[i] == basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) )) ) 
                      ||
                      ( (mu_e[i] == 1) && (mu_k[i] > (1<<mu_j[i])- 1 - basis_->get_numofbw()) ) ))
                    {
                        // nu is not extended right
                        intinfo[i] = 4; // MM
                    }
                    break;
                case 1:
                case 4:
                case 7:
                    break;
                default:
                    abort();
                    break;
            }
        }
        return compute_sum_over_patches (nu_p, intinfo, integralshares);
    }
        
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    double 
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::compute_sum_over_patches (const int nu_p, 
                const MultiIndex<unsigned int, DIM> intinfo, 
                const FixedArray1D<entries,DIM> integralshares) const
    {
        assert (ONEDIMHAARCOUNT == DER_ONEDIMHAARCOUNT); // the case  != is only implemented for the special case DER_COUNT = 1, a=const; you find it in aff_lin_par_eq
        double result =0;
        unsigned int geometry_type;
        int centerpatchnumber;
        MultiIndex<bool, DIM> orientation;
        basis_->get_intersection_geometry(nu_p, intinfo, geometry_type, centerpatchnumber, orientation);
        if (DIM == 2)
        {
            int north, east, northeast;
            if (geometry_type == 3)
            {
                if (centerpatchnumber == -1)
                {
                    north = basis_->get_neighbours(nu_p,0);
                    east = basis_->get_neighbours(nu_p,2);
                    geometry_type = 8;
                }
                else
                {
                    north = basis_->get_neighbours(centerpatchnumber,3);
                    east = basis_->get_neighbours(centerpatchnumber,1);
                    if (north != -1)
                    {
                        northeast = basis_->get_neighbours(north,1);
                        if (northeast == -1)
                        {
                            geometry_type = 10;
                        }
                        else
                        {
                            if (east == -1)
                            {
                                geometry_type = 9;
                            }
                        }
                    }
                    else
                    {
                        northeast = basis_->get_neighbours(east,3);
                        geometry_type = 11;
                    }
                }
            }
            switch (geometry_type)
            {
                case 0: // 1 patch
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a_times_gram, rest;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a_times_gram[eta] = 0;
                        rest[eta] = 0;
                    }
                    if (!orientation[0])
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a_times_gram[y] += agencoeffs_[centerpatchnumber](x,y) * integralshares[0].first[ONEDIMHAARCOUNT-1-x];
                                rest[y] += qgencoeffs_[centerpatchnumber] (x,y) * integralshares[0].first[ONEDIMHAARCOUNT-1-x]
                                        + agencoeffs_[centerpatchnumber] (x,y) * integralshares[0].second[ONEDIMHAARCOUNT-1-x];
                            }
                        }
                    }
                    else
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
    // CLEANUP
                            assert (a_times_gram[y] == 0);
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a_times_gram[y] += agencoeffs_[centerpatchnumber] (x,y) * integralshares[0].first[x];
                                rest[y] += qgencoeffs_[centerpatchnumber] (x,y) * integralshares[0].first[x]
                                        + agencoeffs_[centerpatchnumber] (x,y) * integralshares[0].second[x];
                            }
                        }
                    }
                    if (orientation[1])
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += a_times_gram[y] * integralshares[1].second[y];
                            result += rest[y] * integralshares[1].first[y];
                        }
                    }
                    else
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += a_times_gram[y] * integralshares[1].second[ONEDIMHAARCOUNT-1-y];
                            result += rest[y] * integralshares[1].first[ONEDIMHAARCOUNT-1-y];
                        }
                    }
                    return result;
                }
                    break;
                case 1: // 2 patches: centerpatch and the patch right of it
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a_times_gram, rest;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a_times_gram[eta] = 0;
                        rest[eta] = 0;
                    }
                    const int rightneighbour(basis_->get_neighbours(centerpatchnumber,1));
                    assert (rightneighbour != -1); // assert existence
                    if (!orientation[0])
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[rightneighbour] (x,y)) * integralshares[0].first[x];
                                rest[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[rightneighbour] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[rightneighbour] (x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    else
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[rightneighbour] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[rightneighbour] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[rightneighbour] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    if (orientation[1])
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += a_times_gram[y] * integralshares[1].second[y];
                            result += rest[y] * integralshares[1].first[y];
                        }
                    }
                    else
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += a_times_gram[y] * integralshares[1].second[ONEDIMHAARCOUNT-1-y];
                            result += rest[y] * integralshares[1].first[ONEDIMHAARCOUNT-1-y];
                        }
                    }
                    return result;
                }
                    break;
                case 2: // 2 patches: centerpatch and the one above
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a_times_gram, rest;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a_times_gram[eta] = 0;
                        rest[eta] = 0;
                    }
                    const int upperneighbour(basis_->get_neighbours(centerpatchnumber,3));
                    assert (upperneighbour != -1); // assert existence
                    if (!orientation[1]) // <-- 1 !!
                    {
                        for (int x = 0; x< ONEDIMHAARCOUNT; x++)
                        {
                            for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a_times_gram[x] += (agencoeffs_[centerpatchnumber] (x,ONEDIMHAARCOUNT-1-y) + agencoeffs_[upperneighbour] (x,y)) * integralshares[1].first[y];
                                rest[x] += (qgencoeffs_[centerpatchnumber] (x,ONEDIMHAARCOUNT-1-y) + qgencoeffs_[upperneighbour] (x,y)) * integralshares[1].first[y]
                                        + (agencoeffs_[centerpatchnumber] (x,ONEDIMHAARCOUNT-1-y) + agencoeffs_[upperneighbour] (x,y)) * integralshares[1].second[y];
                            }
                        }
                    }
                    else
                    {
                        for (int x = 0; x< ONEDIMHAARCOUNT; x++)
                        {
                            for (int y = 0 ; y<ONEDIMHAARCOUNT; y++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a_times_gram[x] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[upperneighbour] (x,ONEDIMHAARCOUNT-1-y)) * integralshares[1].first[y];
                                rest[x] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[upperneighbour] (x,ONEDIMHAARCOUNT-1-y)) * integralshares[1].first[y]
                                        + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[upperneighbour] (x,ONEDIMHAARCOUNT-1-y)) * integralshares[1].second[y];
                            }
                        }
                    }
                    if (orientation[0]) // <-- 0 !!
                    {
                        for (unsigned int x=0; x<ONEDIMHAARCOUNT; ++x)
                        {
                            result += a_times_gram[x] * integralshares[0].second[x];
                            result += rest[x] * integralshares[0].first[x];
                        }
                    }
                    else
                    {
                        for (unsigned int x=0; x<ONEDIMHAARCOUNT; ++x)
                        {
                            result += a_times_gram[x] * integralshares[0].second[ONEDIMHAARCOUNT-1-x];
                            result += rest[x] * integralshares[0].first[ONEDIMHAARCOUNT-1-x];
                        }
                    }
                    return result;
                }
                    break;
                case 3: // 4 patches in a square. centerpatch is the lower left one.
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a1_times_gram[eta] = 0;
                        a2_times_gram[eta] = 0;
                        rest1[eta] = 0;
                        rest2[eta] = 0;
                    }
                    
                    if (orientation[0])
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (x,y) + qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    else
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[east] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[northeast] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    if (orientation[1])
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * integralshares[1].second[y];
                            result += (rest1[y]+rest2[ONEDIMHAARCOUNT-1-y]) * integralshares[1].first[y];
                        }
                    }
                    else
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * integralshares[1].second[y];
                            result += (rest1[ONEDIMHAARCOUNT-1-y]+rest2[y]) * integralshares[1].first[y];
                        }
                    }
                    return result;
                }
                    break;
                case 8: // 3 patches, L-shaped support. Southwest is missing, i.e., centerepatch == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a1_times_gram[eta] = 0;
                        a2_times_gram[eta] = 0;
                        rest1[eta] = 0;
                        rest2[eta] = 0;
                    }
                    if (orientation[0])
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (x,y) + qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    else
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[east] (x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[east] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[east] (x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[northeast] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    if (orientation[1])
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * integralshares[1].second[y];
                            result += (rest1[y]+rest2[ONEDIMHAARCOUNT-1-y]) * integralshares[1].first[y];
                        }
                    }
                    else
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * integralshares[1].second[y];
                            result += (rest1[ONEDIMHAARCOUNT-1-y]+rest2[y]) * integralshares[1].first[y];
                        }
                    }
                    return result;
                }
                    break;
                case 9: // 3 patches, L-shaped support. Southeast is missing, i.e., east == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a1_times_gram[eta] = 0;
                        a2_times_gram[eta] = 0;
                        rest1[eta] = 0;
                        rest2[eta] = 0;
                    }
                    
                    if (orientation[0])
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) ) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y) ) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (x,y) ) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (x,y) + qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    else
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[northeast] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    if (orientation[1])
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * integralshares[1].second[y];
                            result += (rest1[y]+rest2[ONEDIMHAARCOUNT-1-y]) * integralshares[1].first[y];
                        }
                    }
                    else
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * integralshares[1].second[y];
                            result += (rest1[ONEDIMHAARCOUNT-1-y]+rest2[y]) * integralshares[1].first[y];
                        }
                    }
                    return result;
                }
                    break;
                case 10: // 3 patches, L-shaped support. Northeast is missing, i.e., northeast == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a1_times_gram[eta] = 0;
                        a2_times_gram[eta] = 0;
                        rest1[eta] = 0;
                        rest2[eta] = 0;
                    }
                    
                    if (orientation[0])
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (x,y) ) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (x,y) ) * integralshares[0].second[x];
                            }
                        }
                    }
                    else
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[east] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].first[x]
                                        + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) ) * integralshares[0].second[x];
                            }
                        }
                    }
                    if (orientation[1])
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * integralshares[1].second[y];
                            result += (rest1[y]+rest2[ONEDIMHAARCOUNT-1-y]) * integralshares[1].first[y];
                        }
                    }
                    else
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * integralshares[1].second[y];
                            result += (rest1[ONEDIMHAARCOUNT-1-y]+rest2[y]) * integralshares[1].first[y];
                        }
                    }
                    return result;
                }
                    break;
                case 11: // 3 patches, L-shaped support. North is missing, i.e., north == -1
                {
                    FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                    for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                    {
                        a1_times_gram[eta] = 0;
                        a2_times_gram[eta] = 0;
                        rest1[eta] = 0;
                        rest2[eta] = 0;
                    }
                    
                    if (orientation[0])
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    else
                    {
                        for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                        {
                            for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                            {
                                // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * integralshares[0].first[x];
                                rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[east] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * integralshares[0].second[x];
                                a2_times_gram[y] += (agencoeffs_[northeast] (x,y)) * integralshares[0].first[x];
                                rest2[y] += (qgencoeffs_[northeast] (x,y)) * integralshares[0].first[x]
                                        + (agencoeffs_[northeast] (x,y)) * integralshares[0].second[x];
                            }
                        }
                    }
                    if (orientation[1])
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * integralshares[1].second[y];
                            result += (rest1[y]+rest2[ONEDIMHAARCOUNT-1-y]) * integralshares[1].first[y];
                        }
                    }
                    else
                    {
                        for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                        {
                            result += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * integralshares[1].second[y];
                            result += (rest1[ONEDIMHAARCOUNT-1-y]+rest2[y]) * integralshares[1].first[y];
                        }
                    }
                    return result;
                }
                    break;
                default:
                    cout << "geometry is not implemented" << endl;
                    abort();
                    break;
            } // end of switch(geometry))
        } // end of DIM == 2
        else
        {
            // DIM == 3. Only Poisson equation implemented!
            if (ONEDIMHAARCOUNT != 1)
            {
                cout << "warning: calling compute_sum_over_patches with DIM >2 and ONEDIMHAARCOUNT > 1.\n This is inefficient. \nThis method computes the Poisson equation, so ONEDIMHAARCOUNT == 1 is sufficient!" << endl;
                abort();
            }
            switch (geometry_type)
            {
                case 0: // 1 patch
                    result = 1;
                    break;
                case 1:
                case 2:
                case 4: // 2 patches
                    result = 2;
                    break;
                case 3:
                case 5:
                case 6: // 4 patches
                    result = 4;
                    break;
                case 7: // 8 patches
                    result = 8;
                    break;
                default:
                    abort();
                    break;
            }
            result *= ( (integralshares[0].second[0] * integralshares[1].first[0]
                        + integralshares[0].first[0] * integralshares[1].second[0])
                        * integralshares[2].first[0]
                      + integralshares[0].first[0] * integralshares[1].first[0]
                      * integralshares[2].second[0]);
            return result;
        }
    }
        
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void 
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::compute_onedim_haar_integrals(const bool reflected,
                const int nui_j,
                const int nui_e,
                const int nui_k,
                const unsigned int nui_basisnum,
                const int mui_j,
                const int mui_e,
                const int mui_k,
                const unsigned int mui_basisnum,
                FixedArray1D<double,ONEDIMHAARCOUNT>& gram,
                FixedArray1D<double,DER_ONEDIMHAARCOUNT>& der) const
    {
        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
        {
            gram[eta] = 0;
        }
        for (unsigned int eta = 0; eta < DER_ONEDIMHAARCOUNT; ++eta)
        {
            der[eta] = 0;
        }
        int nui_k1,nui_k2,mui_k1,mui_k2, nui_der_k1,nui_der_k2,mui_der_k1,mui_der_k2, gram_scale, der_scale;
        basis_->get_bases_infact()[nui_basisnum] -> support(nui_j,nui_e,nui_k,nui_k1,nui_k2);
        basis_->get_bases_infact()[mui_basisnum] -> support(mui_j,mui_e,mui_k,mui_k1,mui_k2);
        unsigned int scale_mui(mui_j+mui_e), scale_nui(nui_j+nui_e);
        unsigned int final_gram_scale(std::max(scale_nui,scale_mui));
        unsigned int final_der_scale(final_gram_scale);
        gram_scale = log2((unsigned int)ONEDIMHAARCOUNT); // assert that wavelet resolution is finer than Haar resolution
        der_scale = log2((unsigned int)DER_ONEDIMHAARCOUNT); // assert that wavelet resolution is finer than Haar resolution
        if ( final_gram_scale < gram_scale)
        {
            final_gram_scale = gram_scale;
        }
        if ( final_der_scale < der_scale)
        {
            final_der_scale = der_scale;
        }

        if (final_der_scale != final_gram_scale)
        {
            nui_der_k1 = nui_k1<<(final_der_scale-scale_nui);
            nui_der_k2 = nui_k2<<(final_der_scale-scale_nui);
            mui_der_k1 = mui_k1<<(final_der_scale-scale_mui);
            mui_der_k2 = mui_k2<<(final_der_scale-scale_mui);
        }
        nui_k1 = nui_k1<<(final_gram_scale-scale_nui);
        nui_k2 = nui_k2<<(final_gram_scale-scale_nui);
        mui_k1 = mui_k1<<(final_gram_scale-scale_mui);
        mui_k2 = mui_k2<<(final_gram_scale-scale_mui);
        
        const unsigned int N_Gauss = basis_->get_bases_infact()[nui_basisnum]->primal_polynomial_degree();
        const double h_gram = ldexp(1.0, -final_gram_scale);
        const double h_der = ldexp(1.0, -final_der_scale);
        Array1D<double> gauss_points, mui_val, mui_der_val, nui_val, nui_der_val;
        double gauss_weight;
        if (reflected)
        {
            int nui_refl_k1, nui_refl_k2;
            
            nui_refl_k1 = (1<<final_gram_scale) - nui_k2;
            nui_refl_k2 = (1<<final_gram_scale) - nui_k1;
            
// CLEANUP        
            assert(nui_refl_k1 < mui_k2); // otherwise no intersection and call to this method useless
            assert(mui_k1 < nui_refl_k2); // otherwise no intersection and call to this method useless
            
            assert(basis_->get_bases_infact()[nui_basisnum]->primal_polynomial_degree() == basis_->get_bases_infact()[mui_basisnum]->primal_polynomial_degree());
            
            if (mui_k1 > nui_refl_k1)
            {
                // mui_k1 is the first patch
                nui_k2 -= (mui_k1-nui_refl_k1); // last relevant patch (w.r.t. the original patch of nui)
            }
            else
            {
                mui_k1 = nui_refl_k1;
            }
            if (mui_k2 < nui_refl_k2)
            {
                // mui_k2 is the last patch
                nui_k1 += (nui_refl_k2-mui_k2); // last relevant patch (w.r.t. the original patch of nui)
            }
            else
            {
                mui_k2=nui_refl_k2;
            }
            
            
            gauss_points.resize(N_Gauss*(mui_k2-mui_k1));
            Array1D<double> gauss_points_refl(N_Gauss*(mui_k2-mui_k1));
            assert ((mui_k2-mui_k1) == (nui_k2-nui_k1));
            // Set up Gauss points and weights for a composite quadrature formula:
            for (int patch = mui_k1, id = 0; patch < mui_k2; patch++) // refers to 2^{-scale}[patch,patch+1]
                for (unsigned int n = 0; n < N_Gauss; n++, id++)
                    gauss_points[id] = h_gram*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
            for (int patch = nui_k1, id = 0; patch < nui_k2; patch++) // refers to 2^{-scale}[patch,patch+1]
                for (unsigned int n = 0; n < N_Gauss; n++, id++)
                    gauss_points_refl[id] = h_gram*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
            
            // - compute point values of the integrands
            if (final_der_scale == final_gram_scale)
            {
                evaluate(*basis_->get_bases_infact()[nui_basisnum], nui_j, nui_e, nui_k, gauss_points_refl, nui_val, nui_der_val);
                evaluate(*basis_->get_bases_infact()[mui_basisnum], mui_j, mui_e, mui_k, gauss_points, mui_val, mui_der_val);
                // - add all integral shares
                if (ONEDIMHAARCOUNT == DER_ONEDIMHAARCOUNT)
                {
                    for (int patch = mui_k1, id1 = 0, id2(nui_val.size()-1); patch < mui_k2; patch++)
                    {
                        // number of the active haar-generator patch, i.e., eta is the biggest integer with
                        // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                        nui_k1 = patch / (1<<(final_gram_scale-gram_scale)); // <-- eta 
                        for (unsigned int n = 0; n < N_Gauss; n++, id1++, id2--) 
                        {
                            gauss_weight = GaussWeights[N_Gauss-1][n] * h_gram;
                            gram[nui_k1] += nui_val[id2] * mui_val[id1] * gauss_weight;
                            // reflection: minus the derivative on nus patch
                            der[nui_k1] -= nui_der_val[id2] * mui_der_val[id1] * gauss_weight;  
                        }
                    }
                }
                else
                {
                    for (int patch = mui_k1, id1_gram = 0, id1_der = 0, id2_gram(nui_val.size()-1), id2_der(nui_der_val.size()-1); patch < mui_k2; patch++)
                    {
                        // number of the active haar-generator patch, i.e., eta is the biggest integer with
                        // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                        nui_k1 = patch / (1<<(final_gram_scale-gram_scale)); // <-- eta 
                        for (unsigned int n = 0; n < N_Gauss; n++, id1_gram++, id2_gram--) 
                        {
                            gauss_weight = GaussWeights[N_Gauss-1][n] * h_gram;
                            gram[nui_k1] += nui_val[id2_gram] * mui_val[id1_gram] * gauss_weight;
                        }
                        nui_der_k1 = patch / (1<<(final_der_scale-der_scale)); // <-- eta 
                        for (unsigned int n = 0; n < N_Gauss; n++, id1_der++, id2_der--) 
                        {
                            gauss_weight = GaussWeights[N_Gauss-1][n] * h_der;
                            // reflection: minus the derivative on nus patch
                            der[nui_der_k1] -= nui_der_val[id2_der] * mui_der_val[id1_der] * gauss_weight;  
                        }
                    }
                }
            }
            else
            {
                evaluate(*basis_->get_bases_infact()[nui_basisnum], 0, nui_j, nui_e, nui_k, gauss_points_refl, nui_val);
                evaluate(*basis_->get_bases_infact()[mui_basisnum], 0, mui_j, mui_e, mui_k, gauss_points, mui_val);
                
                int nui_der_refl_k1 = (1<<final_der_scale) - nui_der_k2;
                int nui_der_refl_k2 = (1<<final_der_scale) - nui_der_k1;
//CLEANUP                
                assert(nui_der_refl_k1 < mui_der_k2); // otherwise no intersection and call to this method useless
                assert(mui_der_k1 < nui_der_refl_k2); // otherwise no intersection and call to this method useless
                
                if (mui_der_k1 > nui_der_refl_k1)
                {
                    // mui_k1 is the first patch
                    nui_der_k2 -= (mui_der_k1-nui_der_refl_k1); // last relevant patch (w.r.t. the original patch of nui)
                }
                else
                {
                    mui_der_k1 = nui_der_refl_k1;
                }
                if (mui_der_k2 < nui_der_refl_k2)
                {
                    // mui_k2 is the last patch
                    nui_der_k1 += (nui_der_refl_k2-mui_der_k2); // last relevant patch (w.r.t. the original patch of nui)
                }
                else
                {
                    mui_der_k2=nui_der_refl_k2;
                }
                
                gauss_points.resize(N_Gauss*(mui_der_k2-mui_der_k1));
                gauss_points_refl.resize(N_Gauss*(mui_der_k2-mui_der_k1));
                assert ((mui_der_k2-mui_der_k1) == (nui_der_k2-nui_der_k1));
                // Set up Gauss points and weights for a composite quadrature formula:
                for (int patch = mui_der_k1, id = 0; patch < mui_der_k2; patch++) // refers to 2^{-scale}[patch,patch+1]
                    for (unsigned int n = 0; n < N_Gauss; n++, id++)
                        gauss_points[id] = h_der*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                for (int patch = nui_der_k1, id = 0; patch < nui_der_k2; patch++) // refers to 2^{-scale}[patch,patch+1]
                    for (unsigned int n = 0; n < N_Gauss; n++, id++)
                        gauss_points_refl[id] = h_der*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                evaluate(*basis_->get_bases_infact()[nui_basisnum], 1, nui_j, nui_e, nui_k, gauss_points_refl, nui_der_val);
                evaluate(*basis_->get_bases_infact()[mui_basisnum], 1, mui_j, mui_e, mui_k, gauss_points, mui_der_val);
                // - add all integral shares
                for (int patch = mui_k1, id1 = 0, id2(nui_val.size()-1); patch < mui_k2; patch++)
                {
                    // number of the active haar-generator patch, i.e., eta is the biggest integer with
                    // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                    nui_k1 = patch / (1<<(final_gram_scale-gram_scale)); // <-- eta 
                    for (unsigned int n = 0; n < N_Gauss; n++, id1++, id2--) 
                    {
                        gauss_weight = GaussWeights[N_Gauss-1][n] * h_gram;
                        gram[nui_k1] += nui_val[id2] * mui_val[id1] * gauss_weight;
                    }
                }
                // - add all integral shares
                for (int patch = mui_der_k1, id1 = 0, id2(nui_der_val.size()-1); patch < mui_der_k2; patch++)
                {
                    // number of the active haar-generator patch, i.e., eta is the biggest integer with
                    // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                    nui_der_k1 = patch / (1<<(final_der_scale-der_scale)); // <-- eta 
                    for (unsigned int n = 0; n < N_Gauss; n++, id1++, id2--) 
                    {
                        gauss_weight = GaussWeights[N_Gauss-1][n] * h_der;
                        // reflection: minus the derivative on nus patch
                        der[nui_der_k1] -= nui_der_val[id2] * mui_der_val[id1] * gauss_weight;  
                    }
                }
            }
            
        }
        else
        {
    // CLEANUP        
            assert(nui_k1 < mui_k2); // otherwise no intersection and call to this method useless
            assert(mui_k1 < nui_k2); // otherwise no intersection and call to this method useless
            assert(basis_->get_bases_infact()[nui_basisnum]->primal_polynomial_degree() == basis_->get_bases_infact()[mui_basisnum]->primal_polynomial_degree());

            mui_k1 = std::max (mui_k1, nui_k1);
            mui_k2 = std::min (mui_k2, nui_k2);
            gauss_points.resize(N_Gauss*(mui_k2-mui_k1));

            // Set up Gauss points and weights for a composite quadrature formula:
            for (int patch = mui_k1, id = 0; patch < mui_k2; patch++) // refers to 2^{-scale}[patch,patch+1]
                for (unsigned int n = 0; n < N_Gauss; n++, id++)
                    gauss_points[id] = h_gram*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;

            if (final_der_scale == final_gram_scale)
            {
                // - compute point values of the integrands
                evaluate(*basis_->get_bases_infact()[nui_basisnum], nui_j, nui_e, nui_k, gauss_points, nui_val, nui_der_val);
                evaluate(*basis_->get_bases_infact()[mui_basisnum], mui_j, mui_e, mui_k, gauss_points, mui_val, mui_der_val);
                // - add all integral shares
                if (ONEDIMHAARCOUNT == DER_ONEDIMHAARCOUNT)
                {
                    for (int patch = mui_k1, id = 0; patch < mui_k2; patch++)
                    {
                        // number of the active haar-generator patch, i.e., eta is the biggest integer with
                        // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                        nui_k1 = patch / (1<<(final_gram_scale-gram_scale));
                        for (unsigned int n = 0; n < N_Gauss; n++, id++) 
                        {
                            gauss_weight = GaussWeights[N_Gauss-1][n] * h_gram;
                            gram[nui_k1] += nui_val[id] * mui_val[id] * gauss_weight;
                            der[nui_k1] += nui_der_val[id] * mui_der_val[id] * gauss_weight;
                        }
                    }
                }
                else
                {
                    for (int patch = mui_k1, id_der = 0, id_gram=0; patch < mui_k2; patch++)
                    {
                        // number of the active haar-generator patch, i.e., eta is the biggest integer with
                        // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                        nui_k1 = patch / (1<<(final_gram_scale-gram_scale));
                        for (unsigned int n = 0; n < N_Gauss; n++, id_gram++) 
                        {
                            gauss_weight = GaussWeights[N_Gauss-1][n] * h_gram;
//                            cout << "gram.size = " << gram.size() << "; nui_val.size() = " << nui_val.soze() << "; mui_val.size() " << mui_val.size() << endl;
                            gram[nui_k1] += nui_val[id_gram] * mui_val[id_gram] * gauss_weight;
                        }
                        nui_der_k1 = patch / (1<<(final_der_scale-der_scale));
                        for (unsigned int n = 0; n < N_Gauss; n++, id_der++) 
                        {
                            gauss_weight = GaussWeights[N_Gauss-1][n] * h_der;
                            der[nui_der_k1] += nui_der_val[id_der] * mui_der_val[id_der] * gauss_weight;
                        }
                    }
                }
            }
            else
            {
                evaluate(*basis_->get_bases_infact()[nui_basisnum], 0, nui_j, nui_e, nui_k, gauss_points, nui_val);
                evaluate(*basis_->get_bases_infact()[mui_basisnum], 0, mui_j, mui_e, mui_k, gauss_points, mui_val);
                // - add all integral shares
                for (int patch = mui_k1, id = 0; patch < mui_k2; patch++)
                {
                    // number of the active haar-generator patch, i.e., eta is the biggest integer with
                    // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                    nui_k1 = patch / (1<<(final_gram_scale-gram_scale));
                    for (unsigned int n = 0; n < N_Gauss; n++, id++) 
                    {
                        gauss_weight = GaussWeights[N_Gauss-1][n] * h_gram;
                        gram[nui_k1] += nui_val[id] * mui_val[id] * gauss_weight;
                    }
                }
                
                mui_der_k1 = std::max (mui_der_k1, nui_der_k1);
                mui_der_k2 = std::min (mui_der_k2, nui_der_k2);
                gauss_points.resize(N_Gauss*(mui_der_k2-mui_der_k1));

                // Set up Gauss points and weights for a composite quadrature formula:
                for (int patch = mui_der_k1, id = 0; patch < mui_der_k2; patch++) // refers to 2^{-scale}[patch,patch+1]
                    for (unsigned int n = 0; n < N_Gauss; n++, id++)
                        gauss_points[id] = h_der*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                evaluate(*basis_->get_bases_infact()[nui_basisnum], 1, nui_j, nui_e, nui_k, gauss_points, nui_der_val);
                evaluate(*basis_->get_bases_infact()[mui_basisnum], 1, mui_j, mui_e, mui_k, gauss_points, mui_der_val);
                // - add all integral shares
                for (int patch = mui_der_k1, id = 0; patch < mui_der_k2; patch++)
                {
                    // number of the active haar-generator patch, i.e., eta is the biggest integer with
                    // 2^{-haar_scale} eta \leq 2^{-scale} k_1
                    nui_der_k1 = patch / (1<<(final_der_scale-der_scale));
                    for (unsigned int n = 0; n < N_Gauss; n++, id++) 
                    {
                        gauss_weight = GaussWeights[N_Gauss-1][n] * h_der;
                        der[nui_der_k1] += nui_der_val[id] * mui_der_val[id] * gauss_weight;
                    }
                }
                
            }
        }
        for (unsigned int eta=0; eta < ONEDIMHAARCOUNT; ++eta)
        {
            if (abs(gram[eta]) < 1e-15)
                gram[eta] = 0;
        }
        for (unsigned int eta=0; eta < DER_ONEDIMHAARCOUNT; ++eta)
        {
            if (abs(der[eta]) < 1e-15)
                der[eta] = 0;
        }
    }
        
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    double
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::norm_A(const unsigned int jmax) const
    {
        return normA_;
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    double
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::norm_Ainv(const unsigned int jmax) const
    {
        return normAinv_;
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::normtest(const unsigned int offset)
    {
        double normA, normAinv;
        std::set<int> Lambda;
        MultiIndex<int, DIM> temp_jmax(basis_->j0()[0]);
        temp_jmax[0]=temp_jmax[0]+offset; 
        unsigned int num(basis_->get_levelnum(temp_jmax, basis_->get_nop()-1));
        Index temp_ind(basis_->last_wavelet(num));
        num = temp_ind.number()+1;
        for (unsigned int n=0;n<num; ++n)
        {
            Lambda.insert(n);
        }
        SparseMatrix<double> A_Lambda;
        setup_stiffness_matrix(*this, Lambda, A_Lambda);
        unsigned int iterations;
        //double help;
        //LanczosIteration(A_Lambda, 1e-6, help, normA, 200, iterations);
        LanczosIteration(A_Lambda, 1e-6, normAinv, normA, 200, iterations);
        normAinv = 1./normAinv;
        cout << "normtest:: offset = " << offset << "; normA = " << normA << "; normAinv = " << normAinv << "; Kondition = " << (normA*normAinv) << endl;
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::RHS(const double eta,                                                
            InfiniteVector<double, Index>& coeffs) const
    {
        coeffs.clear();
        double coarsenorm(0);
        double bound(f_precond_norm_sqr_ - eta*eta);
        //Array1D<std::pair<int,double> > fcoeffs_precond_;
        typename Array1D<std::pair<int, double> >::const_iterator it(fcoeffs_precond_sorted_.begin()), itend(fcoeffs_precond_sorted_.end());
        while ((it != itend) && (coarsenorm < bound))
        {
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(basis_->get_wavelet(it->first), it->second);
            ++it;
        }
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::RHS(const double eta,
            InfiniteVector<double, int>& coeffs) const
    {
        coeffs.clear();
        double coarsenorm(0);
        double bound(f_precond_norm_sqr_ - eta*eta);
        //Array1D<std::pair<int,double> > fcoeffs_precond_;
        typename Array1D<std::pair<int, double> >::const_iterator it(fcoeffs_precond_sorted_.begin()), itend(fcoeffs_precond_sorted_.end());
        while ((it != itend) && (coarsenorm < bound))
        {
            coarsenorm += it->second * it->second;
            coeffs.set_coefficient(it->first, it->second);
            ++it;
        }
    }
      
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::add_ball(const unsigned int& lambdanum,
                                                   Vector<double>& w,
                                                   const int radius,
                                                   const double factor,
                                                   const int maxlevel,
                                                   const CompressionStrategy strategy,
                                                   const bool precond)
    {
        /*
         * For dim=1,2 it is easy to describe all levels in the (1-norm) ball of radius around lambda.
         * (Attention: possibly different minimal levels on different patches!)
         * For higher dimensions this is done recursivly.
         */
        Index temp_lam(basis_->get_wavelet(lambdanum));
        MultiIndex<int,DIM> lam_j(temp_lam.j());
        double d1 = precond ? D(lambdanum) 
                            : 1.0;

        MultiIndex<int,DIM> temp_j;
        unsigned int temp_p;
        if (DIM == 1)
        {
            abort();
            /*
            basis_->get_level(0, temp_j,temp_p);
            //unsigned int start_levelnum(0);
            unsigned int start_num(0);
            if (temp_j[0] < (lam_j[0]-radius))
            {
                temp_j[0] = (lam_j[0]-radius);
                //start_levelnum = basis_->get_levelnum(temp_j,0); // uses: minimal levels differ at most by 1
                start_num = basis_->first_wavelet(basis_->get_levelnum(temp_j,0)).number(); // uses: minimal levels differ at most by 1
                //int start_index = (temp_j[0] >= (lambda.j()[0]-radius))? 0 : basis().first_wavelet(lambda.j()[0]-radius).number();
            }
            temp_j[0] = (lam_j[0]+radius < basis_->get_jmax())? lam_j[0]+radius : basis_->get_jmax();
            //unsigned int end_levelnum(basis_->get_levelnum(temp_j, basis_->get_nop()-1));
            unsigned int end_num(basis_->last_wavelet( basis_->get_levelnum(temp_j, basis_->get_nop()-1)) );
            //int end_index = basis().last_wavelet(min(lambda.j()[0]+radius,basis_->get_jmax())).number();
            double entry;
            for (unsigned int i = start_num; i <= end_num; i++)
            {                
                entry = a(i,lambdanum); // this is too expensive!  But 1d case isn't important
                if (entry != 0)
                {
                    w[i] += precond ? (entry * factor / (d1 * D(i) ))
                                    : entry * factor;
                }
            }
             */
        }
        else if (DIM == 2)
        {
            add_leveldisc_recurse(lambdanum,w,lam_j,-1,radius,factor/d1,list<int>(),false,precond);
#if 0
            
// CLEANUP            
            // cout << "add_ball! DIM = 2" << endl;
            
            // The ball can be described of levellines consisting of levels with the same multidegree
            // The first is determined with the distance of lambda.j and (minimal) j0. The last with basis_->get_jmax()
            // The first level in a levelline is determined with minx = max(j0[0], lambda.j[0]-radius)
            // The last level in a levelline is determined with miny = max(j0[1], lambda.j[1]-radius)
            
            basis_->get_level(0, temp_j,temp_p);
            
            //index_lt j0(this->basis().j0());
            int lambdaline = lam_j[0]+lam_j[1];
            int lowestline = temp_j[0]+temp_j[1];
            int dist2j0=lambdaline-lowestline;
            int dist2maxlevel=basis_->get_jmax()-lambdaline;
            
            MultiIndex<int,DIM> min_j, max_j;
            int xstart,xend,ystart,steps;
            // iterate the levellines. offset relative to lambdas levelline
            unsigned int first_levelnum, last_levelnum; // first and last level on the current levelline
            for (int offset = -std::min(dist2j0,radius); offset < std::min(dist2maxlevel,radius)+1; offset++)
            {
                // iterate over the levels on the levelline

                // ignoring restrictions by j0 for the moment, we have:
                xstart = lam_j[0]-radius+(int)ceil((radius+offset)/2.0); //x coordinate of the first level
                xend = lam_j[0]+(int)floor((radius+offset)/2.0); // same for the last
                ystart = lam_j[1]+(int)floor((radius+offset)/2.0); // and for the second dimension

                // we walk (make steps) through the levelline
                // The first level in the current levelline hast steps=0.
                // restrictions by j0 mean:
                // The first level on the current levelline may have steps >0.
// CLEANUP
                assert (temp_j[0] == basis_->j0()[patch_with_minimal_j0_[0].front()][0]);
                
                bool is_minimal_x(false), is_minimal_y(false);
                if (xstart <= temp_j[0])
                {
                    min_j[0]= temp_j[0]; // xstart+steps;
                    min_j[1]= ystart-(temp_j[0]-xstart);
                    is_minimal_x = true;
                }
                else
                {
                    min_j[0] = xstart;
                    min_j[1] = ystart;
                }
                steps = min(xend-xstart,ystart-basis_->j0()[patch_with_minimal_j0_[1].front()][1]);
                max_j[0]= xstart+steps;
                max_j[1]= ystart-steps;
                
                if (ystart-basis_->j0()[patch_with_minimal_j0_[1].front()][1] <= xend -xstart )
                {
                    is_minimal_y = true;
                }
                
                /*
                 * min_j, max_j contain the minimal and maximal level j on the current levelline
                 * iterate over all those levels. The first level may have the minimal j0[x-direction]. 
                 * In this case we iterate only over those patches with this minimal value.
                 * In the last step the analogous restriction in the y direction applies
                 */
                
// CLEANUP              
//                cout << "offset = " << offset << endl;
//                cout << "lambda = " << temp_lam << "; lambdanum = " << lambdanum << endl;
//                cout << "lambdaline = " << lambdaline << "; ";
//                cout << "lowestline = " << lowestline << "; ";
//                cout << "dist2j0 = " << dist2j0 << "; ";
//                cout << "dist2maxlevel = " << dist2maxlevel << "; ";
//                cout << "basis_->get_jmax() = " << basis_->get_jmax() << "; ";
//                for (unsigned int p=0; p< basis_->get_nop(); ++p)
//                {
//                    cout << "basis_->j0()[" << p << "] = " << basis_->j0()[p] << "; ";
//                }
//                cout << "is_minimal_x = " << is_minimal_x << "; ";
//                cout << "is_minimal_y = " << is_minimal_y << "; ";
//                cout << endl;
//                cout << "patch_with_minimal_j0_[0] = " << endl << "  ";
//                for (list<int>::const_iterator it(patch_with_minimal_j0_[0].begin()), itend(patch_with_minimal_j0_[0].end()); it!= itend; ++it)
//                {
//                    cout << *it << ";  ";
//                }
//                cout << endl;
//                cout << "patch_with_minimal_j0_[1] = " << endl << "  ";
//                for (list<int>::const_iterator it(patch_with_minimal_j0_[1].begin()), itend(patch_with_minimal_j0_[1].end()); it!= itend; ++it)
//                {
//                    cout << *it << ";  ";
//                }
//                cout << endl;
                
                if (is_minimal_x)
                {
                    if (is_minimal_y)
                    {
                        if (steps == 0)
                        {
                            // iterate over patches in the intersection of patch_with_minimal_j0_[0] and patch_with_minimal_j0_[1],
                            // i.e., we are on the lowest possible j0! Only patches where j0()[p] attains the minimal possible value in each dimension
                            // this set might be empty! E.g. j0()[0] = (2,3), j0()[1] = (3,2) => there is no patch with minimal level (2,2)
                            list<int> intersection;
                            set_intersection(patch_with_minimal_j0_[0].begin(), patch_with_minimal_j0_[0].end(),
                                    patch_with_minimal_j0_[1].begin(), patch_with_minimal_j0_[1].end(),
                                    back_inserter(intersection));

        // cleanup
                            for (list<int>::const_iterator it(intersection.begin()), it2(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                            {
                                ++it2;
                                if (it2 != itend)
                                {
                                    assert (*it < *it2);
                                }
                            }
                            assert (min_j == max_j);
                            // it may happen that not all patches between first_levelnum and last_levelnum are in the list "intersection" !!
                            //first_levelnum = basis_->get_levelnum(min_j, intersection.front());
                            //last_levelnum = basis_->get_levelnum(min_j, intersection.back());
                            for (list<int>::const_iterator it(intersection.begin()), it2(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                            {
                                add_level(lambdanum,
                                w,
                                basis_->get_levelnum(min_j,*it),
                                factor/d1,
                                precond);
                            }
                            return;
                        }
                        else
                        {
// TODO : get rid of first_last_levelnum!
sdfsdfsdf                            
                            first_levelnum = basis_->get_levelnum(min_j, patch_with_minimal_j0_[0].front());
                            last_levelnum = basis_->get_levelnum(max_j, patch_with_minimal_j0_[1].back());
                        }
                    }
                    else
                    {
                        // iterate over patches in patch_with_minimal_j0_[0], no restrictions in the y coordinate
                        if (steps == 0)
                        {
sdfsdf                            
                            first_levelnum = basis_->get_levelnum(min_j, patch_with_minimal_j0_[0].front());
                            last_levelnum = basis_->get_levelnum(max_j, patch_with_minimal_j0_[0].back());
                        }
                        else
                        {
sdfsdf                            
                            first_levelnum = basis_->get_levelnum(min_j, patch_with_minimal_j0_[0].front());
                            last_levelnum = basis_->get_levelnum(max_j, basis_->get_nop()-1);
                        }
                    }
                }
                else
                {
                    if (is_minimal_y)
                    {
                        // iterate over the patches in minimal_level_y
                        if (steps == 0)
                        {
sdfsdf                            
                            first_levelnum = basis_->get_levelnum(min_j, patch_with_minimal_j0_[1].front());
                            last_levelnum = basis_->get_levelnum(max_j, patch_with_minimal_j0_[1].back());
                        }
                        else
                        {
sdsdf                            
                            first_levelnum = basis_->get_levelnum(min_j, 0);
                            last_levelnum = basis_->get_levelnum(max_j, patch_with_minimal_j0_[1].back());
                        }
                    }
                    else
                    {
sdfsdf                        
                        // iterate over all patches
                        first_levelnum = basis_->get_levelnum(min_j, 0);
                        last_levelnum = basis_->get_levelnum(max_j, basis_->get_nop()-1);
                    }
                }
                
// CLEANUP
//                cout << "first_levelnum = " << first_levelnum << "; ";
//                cout << "last_levelnum = " << last_levelnum << endl;
                
                for (unsigned int levelnum = first_levelnum; levelnum <= last_levelnum; levelnum++)
                {
                    //cout << "adding levelnum = " << levelnum << endl;
//                    cout << "cached_qtproblem.cpp :: levelnum = " << levelnum << "; pre add_level: w[" << M_ << "] = " << w[M_] << endl;
//                    cout << "cached_qtproblem.cpp :: factor = " << factor << "; d1 = " << d1 << "; factor/d1 = " <<  factor/d1 << endl;
                    add_level(lambdanum,
                            w,
                            levelnum,
                            factor/d1,
                            precond);
                    // not very efficient. different posibility: precompute in compose_wavelets (then an additional parameter "precond" needs to be added there)
//                    if  (precond)
//                    {
//                        for (unsigned int n=basis_->first_wavelet(levelnum).number(); n<=basis_->last_wavelet(levelnum).number(); ++n)
//                        {
//// CLEANUP
////                            if (n == M_)
////                            {
////                                cout << "cached_qtproblem.cpp :: w[" << n << "] = " << w[n] << "; D(" << n << ") = " << D(n) << endl;
////                            }
//                            if (w[n] != 0)
//                            {
//                                w[n] = w[n] / D(n) ;
//                            }
//                            if (abs (w[n] - factor*a(n,lambdanum)/D(lambdanum)/D(n)) > 1e-14 )
//                            {
//                                cout << "lambdanum = " << lambdanum << "; lambda = " << *basis_->get_wavelet(lambdanum) << "; n = " << n << "; mu = " << *basis_->get_wavelet(n) << "; w[" << n << "] = " << w[n] << "; f*a(" << n << ", " << lambdanum << ") /D(" << lambdanum << ")/D(" << n << ")= " << factor*a(n,lambdanum)/D(lambdanum)/D(n) << endl;
//                                cout << "diff = " << (abs(factor*a(n,lambdanum)/D(lambdanum)/D(n) - w[n])) << "; quot = " << (factor*a(n,lambdanum)/D(lambdanum)/D(n) / w[n]) << endl;
//                            }
//                            assert (abs (w[n] - factor*a(n,lambdanum)/D(lambdanum)/D(n)) < 1e-14 );
//                        }
//                    }
//                    else
//                    {
//                        for (unsigned int n=basis_->first_wavelet(levelnum).number(); n<=basis_->last_wavelet(levelnum).number(); ++n)
//                        {
//                            if (w[n] != factor*a(n,lambdanum))
//                            {
//                                cout << "lambdanum = " << lambdanum << "; w[" << n << "] = " << w[n] << "; f*a(" << n << ", " << lambdanum << ") = " << factor*a(n,lambdanum) << endl;
//                                cout << "diff = " << (abs(factor*a(n,lambdanum) - w[n])) << "; quot = " << (factor*a(n,lambdanum) / w[n]) << endl;
//                            }
//                            assert (w[n] == factor*a(n,lambdanum));
//                        }
//                    }
////                    cout << "cached_qtproblem.cpp :: post add_level: w[" << M_ << "] = " << w[M_] << endl;
                }
                
            }
#endif              
        }
        else // dim > 2. iteration over all levels in 'range' is done recursivly for arbitrary dimensions
        {
            // iterate over all possible values for the first component of the level j
            // for each j[0]: add levels with valid values for j[1],...,j[DIM-1]
            // this can be done recursively
            // for the last two dimensions instead of doing a final recursion, proceed as for the case DIM=2 (above in this method)
            
            int xstart, xend;
            xstart = lam_j[0]-radius;
            bool set_is_restriction(false);
            // basis_->j0()[patch_with_minimal_j0_[0].first][0] == lowest possible value for j[0]
            if (xstart <= basis_->j0()[patch_with_minimal_j0_[0].front()][0])
            {
                xstart = basis_->j0()[patch_with_minimal_j0_[0].front()][0];
                set_is_restriction =(patch_with_minimal_j0_[0].size() != basis_->get_nop()); // check whether all patches allow for the lowest possible value for j[0]
            }
            temp_j = lam_j;
            assert (radius >= 0);
            assert (basis_->get_jmax() - multi_degree(lam_j) >= 0);
            xend = lam_j[0] + std::min (radius, (int) (basis_->get_jmax() - multi_degree(lam_j)));
            
            temp_j[0] = xstart;
            add_leveldisc_recurse(lambdanum,w,temp_j,0,radius-abs(xstart-lam_j[0]),factor/d1,patch_with_minimal_j0_[0],set_is_restriction,precond);
                
            for (unsigned int x = xstart+1; x<= xend; ++x)
            {
                // first component of the current level is set to "x"
                // add all possible levels with this value
                // ASSUMPTION: minimal levels differ at most by 1 => for x>xstart the value x is possible for each patch
                temp_j[0] = x;
                add_leveldisc_recurse(lambdanum,w,temp_j,0,radius-abs(x-lam_j[0]),factor/d1,list<int>(),false,precond);
            }
            
        }
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::add_leveldisc_recurse(
                const unsigned int& lambdanum,
                Vector<double>& w,
                const MultiIndex<int,DIM> center_j,
                const int fixed_dimensions,
                const unsigned int radius,
                const double factor,
                const list<int> & j0_restriction_of_patches,
                const bool set_is_restriction,
                const bool precond)
    {
        assert (DIM-fixed_dimensions>2);
        MultiIndex<int,DIM> temp_j;
        int xstart, xend;
        list<int> intersection;
        
        if (DIM-fixed_dimensions > 3)
        {
            xstart = center_j[fixed_dimensions+1]-radius;
            bool restricted(set_is_restriction);
            // basis_->j0()[patch_with_minimal_j0_[0].first][0] == lowest possible value for j[0]
            if (xstart <= basis_->j0()[patch_with_minimal_j0_[fixed_dimensions+1].front()][fixed_dimensions+1])
            {
                xstart = basis_->j0()[patch_with_minimal_j0_[fixed_dimensions+1].front()][fixed_dimensions+1];
                if (patch_with_minimal_j0_[fixed_dimensions+1].size() != basis_->get_nop())
                {
                    if (restricted)
                    {
                        set_intersection(patch_with_minimal_j0_[fixed_dimensions+1].begin(), patch_with_minimal_j0_[fixed_dimensions+1].end(),
                                j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                back_inserter(intersection));
                        if (intersection.empty())
                        {
                            return;
                        }
                    }
                    else
                    {
                        restricted = true;
                        intersection = patch_with_minimal_j0_[fixed_dimensions+1];
                    }
                }
            }
            temp_j = center_j;
            assert (radius >= 0);
            assert (basis_->get_jmax() - multi_degree(center_j) >= 0);
            xend = center_j[fixed_dimensions+1] + std::min ((int) radius, (int) (basis_->get_jmax()) - (int)multi_degree(center_j));
            
            temp_j[fixed_dimensions+1] = xstart;
            add_leveldisc_recurse(lambdanum,w,temp_j,fixed_dimensions+1,radius-abs(xstart-center_j[fixed_dimensions+1]),factor,  intersection,restricted,precond);
            for (unsigned int x = xstart+1; x<= xend; ++x)
            {
                // first component of the current level is set to "x"
                // add all possible levels with this value
                // ASSUMPTION: minimal levels differ at most by 1 => for x>xstart the value x is possible for each patch
                temp_j[fixed_dimensions+1] = x;
                add_leveldisc_recurse(lambdanum,w,temp_j,fixed_dimensions+1,radius-abs(x-center_j[fixed_dimensions+1]),factor,  j0_restriction_of_patches,set_is_restriction,precond);
            }
        }
        else
        {
            // there are only 2 dimensions left
            // if set_is_restriction == false: the code to iterate over the last 2 dimensions looks the same as in add_ball for DIM==2
            // else: the current levelline has holes, i.e., only patches in the restriction set may be active
            
            // The ball can be described of levellines consisting of levels with the same multidegree
            // The first is determined with the distance of lambda.j and (minimal) j0. The last with basis_->get_jmax()
            // The first level in a levelline is determined with minx = max(j0[DIM-2], lambda.j[DIM-2]-radius)
            // The last level in a levelline is determined with miny = max(j0[DIM-1], lambda.j[DIM-1]-radius)
            
            MultiIndex<int, 2> min_j, max_j;
            unsigned int temp_p;
            basis_->get_level(0, temp_j,temp_p);
            
            //index_lt j0(this->basis().j0());
            int lambdaline = multi_degree(center_j);
            int lowestline = multi_degree(temp_j);
            int dist2j0=lambdaline-lowestline;
            int dist2maxlevel=basis_->get_jmax()-lambdaline;
            int ystart,temp_i;
            bool is_minimal_x(false), is_minimal_y(false); // is_minimal_x == true means that the first level of the current levelline (in the current ball) hits j0[DIM-2]
            list<int> intersection_x, intersection_y; // intersection_x contains (if needed) the intersection of j0_restriction_of_patches and patch_with_minimal_j0_[DIM-2]  
            
            //bool j0_restricts_x(set_is_restriction), j0_restrics_y(set_is_restriction); // false means that all patches yield valid levels
            //List<int>* list_x; // if j0_restricts_x == true this points to the list with the valid patches
            //list<int>* list_y;
            //unsigned int first_levelnum, last_levelnum; // first and last level on the current levelline
            
            
            // iterate the levellines. offset relative to center_j's levelline
            for (int offset = -std::min(dist2j0,(int)radius); offset < std::min(dist2maxlevel,(int)radius)+1; offset++)
            {
                // iterate over the levels on the levelline
                // ignoring restrictions by j0 for the moment, we have:
                
                xstart = center_j[DIM-2]-radius+(int)ceil((radius+offset)/2.0); //x coordinate of the first level
                xend = center_j[DIM-2]+(int)floor((radius+offset)/2.0); // same for the last
                ystart = center_j[DIM-1]+(int)floor((radius+offset)/2.0); // and for the second dimension
                // we walk (make steps) through the levelline
                // The level (xstart,ystart) in the current levelline hast steps=0.
                // restrictions by j0[DIM-2] or j0[DIM-1] mean:
                // The first level on the current levelline may have steps >0.
                
                // ignoring the restriction by set_is_restriction for the moment, we have
                
                temp_i = basis_->j0()[patch_with_minimal_j0_[DIM-2].front()][DIM-2];
                if (xstart <= temp_i)
                {
                    min_j[0]= temp_i;
                    min_j[1]= ystart-(temp_i-xstart);
                    is_minimal_x = true; // restrictions by patch_with_minimal_j0_[DIM-2] apply
                }
                else
                {
                    min_j[0] = xstart;
                    min_j[1] = ystart;
                }
                temp_i = min(xend-xstart,ystart-basis_->j0()[patch_with_minimal_j0_[DIM-1].front()][DIM-1]);
                max_j[0]= xstart+temp_i;
                max_j[1]= ystart-temp_i;

                if (ystart-basis_->j0()[patch_with_minimal_j0_[DIM-1].front()][DIM-1] <= xend - xstart )
                {
                    is_minimal_y = true; // restrictions by patch_with_minimal_j0_[DIM-1] apply
                }
                
                // how many j's are legal between xmin and xmax?
                temp_i = max_j[0] - min_j[0];
                if (temp_i < 0)
                {
                    return;
                }
// cleanup
                assert (temp_i == (min_j[1]-max_j[1]));
                // incorporate set_is_restriction
                temp_j = center_j;
                temp_j[DIM-2] = min_j[0];
                temp_j[DIM-1] = min_j[1];          

                if (set_is_restriction)
                {
                    // restrictions by j0_restriction_of_patches apply
                    if (is_minimal_x)
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] apply
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(patch_with_minimal_j0_[DIM-2].begin(), patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
                                set_intersection(patch_with_minimal_j0_[DIM-1].begin(), patch_with_minimal_j0_[DIM-1].end(),
                                        intersection_x.begin(), intersection_x.end(),
                                        back_inserter(intersection));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection.begin()), it2(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection and with j = min_j
                                for (list<int>::const_iterator it(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                set_intersection(patch_with_minimal_j0_[DIM-2].begin(), patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
                                set_intersection(patch_with_minimal_j0_[DIM-1].begin(), patch_with_minimal_j0_[DIM-1].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_y));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_x.begin()), it2(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                for (list<int>::const_iterator it(intersection_y.begin()), it2(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                // iterate over all levels with p\in intersection_x and with j = min_j
                                for (list<int>::const_iterator it(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                for (int i=1; i<temp_i;++i)
                                {
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                    for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            basis_->get_levelnum(temp_j,*it),
                                            factor,
                                            precond);
                                    }
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(patch_with_minimal_j0_[DIM-2].begin(), patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_x.begin()), it2(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection_x.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection and with j = min_j
                                for (list<int>::const_iterator it(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                set_intersection(patch_with_minimal_j0_[DIM-2].begin(), patch_with_minimal_j0_[DIM-2].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_x));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_x.begin()), it2(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                // iterate over all levels with p\in intersection_x and with j = min_j
                                for (list<int>::const_iterator it(intersection_x.begin()), itend(intersection_x.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                for (int i=1; i<=temp_i;++i)
                                {
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                    for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            basis_->get_levelnum(temp_j,*it),
                                            factor,
                                            precond);
                                    }
                                }
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                            }
                        }
                    }
                    else
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] do not apply                        
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(patch_with_minimal_j0_[DIM-1].begin(), patch_with_minimal_j0_[DIM-1].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_y));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_y.begin()), it2(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection_y.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection_y and with j = min_j
                                for (list<int>::const_iterator it(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                set_intersection(patch_with_minimal_j0_[DIM-1].begin(), patch_with_minimal_j0_[DIM-1].end(),
                                        j0_restriction_of_patches.begin(), j0_restriction_of_patches.end(),
                                        back_inserter(intersection_y));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection_y.begin()), it2(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                // iterate over all levels with p\in j0_restriction_of_patches and with j = min_j, steps = 0,...,temp_i-1
                                for (int i=0; i<temp_i;++i)
                                {
                                    for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            basis_->get_levelnum(temp_j,*it),
                                            factor,
                                            precond);
                                    }
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                }
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(intersection_y.begin()), itend(intersection_y.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply

                            // iterate over all levels with p\in intersection_x and with j = min_j
                            for (int i=0; i<=temp_i;++i)
                            {
                                for (list<int>::const_iterator it(j0_restriction_of_patches.begin()), itend(j0_restriction_of_patches.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
                            }
// CLEANUP
                            assert ((temp_j[DIM-2]-1) == max_j[0]);
                            assert ((temp_j[DIM-1]+1) == max_j[1]);
                        }
                    }
                }
                else
                {
                    // restrictions by j0_restriction_of_patches do not apply
                    if (is_minimal_x)
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] apply
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                set_intersection(patch_with_minimal_j0_[DIM-2].begin(), patch_with_minimal_j0_[DIM-2].end(),
                                        patch_with_minimal_j0_[DIM-1].begin(), patch_with_minimal_j0_[DIM-1].end(),
                                        back_inserter(intersection));
// CLEANUP assert sorting
                                for (list<int>::const_iterator it(intersection.begin()), it2(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    ++it2;
                                    if (it2 != itend)
                                    {
                                        assert (*it < *it2);
                                    }
                                }
                                if (intersection.empty())
                                {
                                    return;
                                }
                                // iterate over all levels with p\in intersection and with j = min_j
                                for (list<int>::const_iterator it(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
//                                cout << "start: temp_j = " << temp_j << endl;
                                for (list<int>::const_iterator it(patch_with_minimal_j0_[DIM-2].begin()), itend(patch_with_minimal_j0_[DIM-2].end()); it!= itend; ++it)
                                {
//                                    cout << "start: adding level temp_j = " << temp_j << "; p = " << *it << endl;
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                                for (int i=1; i<temp_i;++i)
                                {
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
//                                    cout << "center: temp_j = " << temp_j << endl;
                                    for (unsigned int p=0; p<basis_->get_nop();++p)
                                    {
//                                        cout << "center: adding level temp_j = " << temp_j << "; p = " << p << endl;
                                        add_level(lambdanum,
                                            w,
                                            basis_->get_levelnum(temp_j,p),
                                            factor,
                                            precond);
                                    }
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
//                                cout << "end: temp_j = " << temp_j << endl;
//                                cout << "min_j = " << min_j << endl;
//                                cout << "max_j = " << max_j << endl;
                                
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(patch_with_minimal_j0_[DIM-1].begin()), itend(patch_with_minimal_j0_[DIM-1].end()); it!= itend; ++it)
                                {
//                                    cout << "end: adding level temp_j = " << temp_j << "; p = " << *it << endl;
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            } 
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply
                            // iterate over all levels with p\in intersection_x and with j = min_j
                            for (list<int>::const_iterator it(patch_with_minimal_j0_[DIM-2].begin()), itend(patch_with_minimal_j0_[DIM-2].end()); it!= itend; ++it)
                            {
                                add_level(lambdanum,
                                    w,
                                    basis_->get_levelnum(temp_j,*it),
                                    factor,
                                    precond);
                            }
                            for (int i=1; i<=temp_i;++i)
                            {
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
                                for (unsigned int p=0; p<basis_->get_nop();++p)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,p),
                                        factor,
                                        precond);
                                }
                            }
// CLEANUP
                            assert (temp_j[DIM-2] == max_j[0]);
                            assert (temp_j[DIM-1] == max_j[1]);
                        }
                    }
                    else
                    {
                        // restrictions by patch_with_minimal_j0_[DIM-2] do not apply
                        if (is_minimal_y)
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] apply
                            if (temp_i == 0)
                            {
                                // all levels in the levelline have the same j
                                // iterate over all levels with p\in intersection_y and with j = min_j
                                for (list<int>::const_iterator it(patch_with_minimal_j0_[DIM-1].begin()), itend(patch_with_minimal_j0_[DIM-1].end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                            else
                            {
                                // min_j != max_j
                                // iterate over all levels with p\in j0_restriction_of_patches and with j = min_j, steps = 0,...,temp_i-1
                                for (int i=0; i<temp_i;++i)
                                {
                                    for (unsigned int p=0; p<basis_->get_nop();++p)
                                    {
                                        add_level(lambdanum,
                                            w,
                                            basis_->get_levelnum(temp_j,p),
                                            factor,
                                            precond);
                                    }
                                    ++temp_j[DIM-2];
                                    --temp_j[DIM-1];
                                }
// CLEANUP
                                assert (temp_j[DIM-2] == max_j[0]);
                                assert (temp_j[DIM-1] == max_j[1]);
                                for (list<int>::const_iterator it(patch_with_minimal_j0_[DIM-1].begin()), itend(patch_with_minimal_j0_[DIM-1].end()); it!= itend; ++it)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,*it),
                                        factor,
                                        precond);
                                }
                            }
                        }
                        else
                        {
                            // restrictions by patch_with_minimal_j0_[DIM-1] do not apply

                            // iterate over all levels with p\in intersection_x and with j = min_j
                            for (int i=0; i<=temp_i;++i)
                            {
                                for (unsigned int p=0; p<basis_->get_nop();++p)
                                {
                                    add_level(lambdanum,
                                        w,
                                        basis_->get_levelnum(temp_j,p),
                                        factor,
                                        precond);
                                }
                                ++temp_j[DIM-2];
                                --temp_j[DIM-1];
                            }
// CLEANUP
                            assert ((temp_j[DIM-2]-1) == max_j[0]);
                            assert ((temp_j[DIM-1]+1) == max_j[1]);
                        }
                    }
                }
            }
        }

}
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::add_level(const unsigned int& lambdanum,
                Vector<double>& w,
                const unsigned int levelnum,
                const double factor,
                const bool precond)
    {
        // for fixed lambda: compute entries a(mu_2,lambda) for all mu_2 on level levelnum that intersect lambda
        // part 1: for each dimension: extract 1d Block a(mu2_i, lam_i) of 1d integrals with |mu2_i| = |levelnum_j[i]|
        // as in a() these blocks may or may not already be in the cache. So maybe they have to be computed here
        // part 2: compose all relevant entries of w from the extracted 1d Blocks

        // Part 1: extract 1d Blocks from the cache (similar to a() ))
        Index lambda (basis_->get_wavelet(lambdanum));
        const MultiIndex<int, DIM> lambda_j(lambda.j()), lambda_e(lambda.e()), lambda_k(lambda.k());
        const unsigned int lambda_p(lambda.p());
        MultiIndex<unsigned int, DIM> intinfo;
        MultiIndex<int, DIM> mu_j;
        unsigned int mu_p;
        basis_->get_level(levelnum,mu_j,mu_p);
        //bool temp_b(basis_->get_LMR_info(lambdanum,mu_j,mu_p, intinfo));
        bool temp_b(basis_->get_LMR_info(lambdanum,mu_p, intinfo));
        if (!temp_b) 
            return;
        FixedArray1D<int,DIM> kmingen, kmaxgen, kminwav, kmaxwav;
        unsigned int lami_basisnum, mui_basisnum;
        FixedArray1D<bool,DIM> mu_min_type;
        FixedArray1D< Block*, DIM> waveletBlock;
        FixedArray1D< Block*, DIM> generatorBlock;
        for (unsigned int i=0; i<DIM; i++)
        {
            lami_basisnum = (((basis_->get_bc()[lambda_p][2*i])?0:2) + ((basis_->get_bc()[lambda_p][2*i+1])?0:1));
            mui_basisnum = (((basis_->get_bc()[mu_p][2*i])?0:2) + ((basis_->get_bc()[mu_p][2*i+1])?0:1));
            mu_min_type[i] = (mu_j[i]==basis_->j0()[mu_p][i]) ? true:false;
            ColumnCache* cachePointer;
            Column* cachePointerGen;
            switch (intinfo[i])
            {
                case 0:
                case 4:
                case 8: // LL MM RR -> type 1
                    cachePointer = & typeIcache_[4*lami_basisnum + mui_basisnum];
                    if (mu_min_type[i])
                        cachePointerGen = & typeIcachedGens_[4*lami_basisnum + mui_basisnum];
                    break;
                case 1:
                case 7: 
                case 2:
                case 6: // LM RM LR RL -> type 2
                    cachePointer = &typeIIcache_[4*(lami_basisnum-1) + mui_basisnum];
                    if (mu_min_type[i])
                        cachePointerGen = & typeIIcachedGens_[4*(lami_basisnum-1) + mui_basisnum];
                    break;
                case 3: // ML -> type 3
                    // EXTREME CAUTION:: mui_basisnum - (intinfo[i] == 5)?0:1 produces garbage!
                    assert (mui_basisnum != 1);
                    cachePointer = &typeIIIcache_[4*lami_basisnum + mui_basisnum - 1];
                    if (mu_min_type[i])
                        cachePointerGen = & typeIIIcachedGens_[4*lami_basisnum + mui_basisnum - 1];
                    break;
                case 5: // MR -> type 3
                    assert (mui_basisnum != 2);
                    cachePointer = &typeIIIcache_[4*lami_basisnum + ((mui_basisnum == 1)?0:3)];
                    if (mu_min_type[i])
                        cachePointerGen = & typeIIIcachedGens_[4*lami_basisnum + ((mui_basisnum == 1)?0:3)];
                    break;
                default:
                    abort();
            }
            // search for column 'lami'
            typename QTBASIS::IntervalBasis::Index lami(lambda_j[i],lambda_e[i],lambda_k[i],basis_->get_bases_infact()[lami_basisnum]);
            unsigned int lami_num = lami.number();
            typename ColumnCache::iterator col_lb(cachePointer->lower_bound(lami_num));
            typename ColumnCache::iterator col_it(col_lb);
            if (col_lb == cachePointer->end() ||
                cachePointer->key_comp()(lami_num, col_lb->first))
            {
                // the lami-th column has not been requested so far
                // insert a new column and continue with the blocks
                typedef typename ColumnCache::value_type value_type;
                col_it = cachePointer->insert(col_lb, value_type(lami_num, Column()));
            }
            // col_it points to the column of lami
            Column& col(col_it->second);
            // check whether the level 'mu_i' belongs to has already been calculated
            int mui_levelnum (mu_j[i] - basis_->j0()[mu_p][i]);
            typename Column::iterator lb(col.lower_bound(mui_levelnum));
            typename Column::iterator it(lb);
            if (lb == col.end() ||
                col.key_comp()(mui_levelnum, lb->first))
            {
                // no entries have ever been computed for this column and this level
                // insert a new block, 
                // and then compute the whole level block 
                //      (\int psi_mui psi_lami)_mui, (\int psi_mui' psi_lami')_mui
                // for all mui that intersect lami. 
                // if mui is on the minimal level we also need to add a new Block() to
                // the generator integral cache, i.e., we need to 
                // integrate against all generators on the lowest level
                typedef typename Column::value_type value_type;
                it = col.insert(lb, value_type(mui_levelnum, Block()));
                waveletBlock[i] = & it->second;
                bool gen_intersection_i, wav_intersection_i;
                basis_->get_onedim_intersections(intinfo[i],
                        lambda_j[i],
                        lambda_e[i],
                        lambda_k[i],
                        lami_basisnum,
                        (mui_levelnum == 0),
                        mu_j[i],
                        mui_basisnum,
                        kmingen[i],
                        kmaxgen[i],
                        kminwav[i],
                        kmaxwav[i],
                        gen_intersection_i,
                        wav_intersection_i);
                FixedArray1D<double,ONEDIMHAARCOUNT> gram;
                FixedArray1D<double,DER_ONEDIMHAARCOUNT> der;
                if (mu_min_type[i])
                {
                    typename Column::iterator lb2(cachePointerGen->lower_bound(lami_num));
                    typename Column::iterator it2(lb2);
                    if (lb2 == cachePointerGen->end() ||
                        cachePointerGen->key_comp()(lami_num, lb2->first))
                    {
                        // no entries have ever been computed for this column and this level
                        // insert a new block, 
                        // and then compute the whole level block 
                        //      (\int psi_nui psi_mui)_nui, (\int psi_nui' psi_mui')_nui
                        // for all nui that intersect mui. 
                        // if mui is on the minimal level we also need to add a new Block() to
                        // the generator integral cache, i.e., we need to 
                        // integrate against all generators on the lowest level
                        typedef typename Column::value_type value_type;
                        it2 = cachePointerGen->insert(lb2, value_type(lami_num, Block()));
                        generatorBlock[i] =  & it2->second;
                        //Block& block2(it2->second);
                        if (gen_intersection_i)
                        {
                            for (int kgen = kmingen[i]; kgen <= kmaxgen[i]; ++kgen)
                            {
                                compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), lambda_j[i], lambda_e[i], lambda_k[i], lami_basisnum,
                                        mu_j[i], 0, kgen, mui_basisnum,
                                        gram,
                                        der);
                                typedef typename Block::value_type value_type_block;
                                generatorBlock[i]->insert(generatorBlock[i]->end(), value_type_block(kgen, make_pair(gram,der)));
                            }
                        }
                    }
                    else
                    {
                        cout << "just inserted a wavelet block, but gen block was already existent?!!" << endl;
                        abort();
                    }
                }
                // code invariant: there exists a Block() (maybe empty) corresponding to lami and mui
                
                // if (no intersections) => current block (wavCache) remains empty
                // if we are on the lowest possible level for mui, we add an empty Block to the genCache
                // advantage: genCache[laminum] exists and contains a block. 
                // code invariant: to access an entry: load the block and check if there is an entry (simple! same whether there is an entry or not!)
                // it would be possible to simply leave the genCache as it is.
                // advantage: nothing to do at this point
                // access to an entry: check whether there is a column in genCache. If not: value 0. If there is: get Block and take value from it.
                // This leads to 2 checks per call to a() instead of one. 
                // This is cheaper than method 1 above if lambda does not intersect at all with basis functions mu.
                // However, most lambdas will have some intersection? In this case this variant leads to 1 additional check
                // so ... hopefully this makes the average access time to the cache faster
                
                if (!wav_intersection_i) 
                    return;
                // wav_intersection_i == true guarantees kminwavi <=kmaxwavi and that the values are meaningful
                for (int kwav = kminwav[i]; kwav <= kmaxwav[i]; ++kwav)
                {
                    compute_onedim_haar_integrals(!(((intinfo[i] == 0) || (intinfo[i] == 4) ) || (intinfo[i] == 8) ), lambda_j[i], lambda_e[i], lambda_k[i], lami_basisnum,
                            mu_j[i], 1, kwav, mui_basisnum,
                            gram,
                            der);

                    typedef typename Block::value_type value_type_block;
                    waveletBlock[i]->insert(waveletBlock[i]->end(), value_type_block(kwav, make_pair (gram,der) ));
                }
            }
            else // column and levelblock already exist
            {
                waveletBlock[i] = & it->second;
                if (waveletBlock[i]->size() == 0)
                    return;
                if (mu_min_type[i])
                {
                    typename Column::iterator lb2(cachePointerGen->lower_bound(lami_num));
                    generatorBlock[i] = & lb2->second;
                }
            }
        } // end of loop over dim
        // part 2: compose all relevant entries of w from the extracted 1d Blocks
        // part 2a: include information about mu into LMR info:
        FixedArray1D<Array1D<unsigned int> ,DIM> mu_gen_adapted_intinfo;
        FixedArray1D<Array1D<unsigned int> ,DIM> mu_wav_adapted_intinfo;
        for (unsigned int i=0; i<DIM; ++i)
        {
            kminwav[i] = waveletBlock[i]->begin()->first;
            kmaxwav[i] = waveletBlock[i]->rbegin()->first;
            assert (waveletBlock[i]->size() == (kmaxwav[i]-kminwav[i]+1));
            if (mu_min_type[i] && (generatorBlock[i]->size() != 0))
            {
                // loop over all intersecting generators in current dimension
                mu_gen_adapted_intinfo[i].resize(generatorBlock[i]->size());
                kmingen[i] = generatorBlock[i]->begin()->first;
                kmaxgen[i] = generatorBlock[i]->rbegin()->first;
                assert (generatorBlock[i]->size() == (kmaxgen[i]-kmingen[i]+1));
                typename Block::iterator block_lb(generatorBlock[i]->lower_bound(kmingen[i]));
                if (block_lb == generatorBlock[i]->end() ||
                        generatorBlock[i]->key_comp()(kmingen[i], block_lb->first))
                {
                    cout << "error! 1d cache contains generator block for current level but valid entry kmin_gen[i] is missing!" << endl;
                    abort();
                }
                for (int k = kmingen[i], n(0); k<= kmaxgen[i]; ++k, ++n)
                {
                    // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
                    mu_gen_adapted_intinfo[i][n] = intinfo[i];
                    switch (intinfo[i])
                    {
                        case 0:
                            // LL: check whether mu is a left boundary gen/wav
                            if (!(k == basis_-> get_bases(mu_p,i)->DeltaLmin()) )
                            {
                                // nu is not extended left
                                mu_gen_adapted_intinfo[i][n] = 4; // MM
                            }
                            break;
                        case 2:
                            // LR: check whether mu is a right boundary gen/wav
                            if (!(k == basis_-> get_bases(mu_p,i)->DeltaRmax(  basis_->get_j0(mu_p,i) )) )
                              
                            {
                                mu_gen_adapted_intinfo[i][n] = 1; // LM
                            }
                            break;
                        case 3:
                            // ML : check whether mu is a left boundary gen/wav
                            if (!(k == basis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                            {
                                // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                                abort();
                            }
                            break;
                        case 5:
                            // MR: check whether mu is a right boundary gen/wav
                            if (!(k == basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) )) )
                            {
                                // nu is not extended right. There is no intersection of mu and nu! This should not happen at this point
                                abort();
                            }
                            break;
                        case 6:
                            // RL : check whether mu is a left boundary gen/wav
                            if (!(k == basis_-> get_bases(mu_p,i)->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                            {
                                // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                                mu_gen_adapted_intinfo[i][n] = 7; // RM
                            }
                            break;
                        case 8:
                            // RR: check whether mu is a right boundary gen/wav
                            if (!(k == basis_-> get_bases(mu_p,i)->DeltaRmax( basis_->get_j0(mu_p,i) )) ) 
                            {
                                // nu is not extended right
                                mu_gen_adapted_intinfo[i][n] = 4; // MM
                            }
                            break;
                        case 1:
                        case 4:
                        case 7:
                            break;
                        default:
                            abort();
                            break;
                    }
                } // end of loop over generatorBlock
            }
            // loop over all intersecting wavelets in current dimension
            mu_wav_adapted_intinfo[i].resize(kmaxwav[i]-kminwav[i]+1);
            typename Block::iterator block_lb(waveletBlock[i]->lower_bound(kminwav[i]));
            if (block_lb == waveletBlock[i]->end() ||
                    waveletBlock[i]->key_comp()(kminwav[i], block_lb->first))
            {
                cout << "error! 1d cache contains wavelet block for current level but valid entry kmin_wav[i] is missing!" << endl;
                abort();
            }
            for ( int k = kminwav[i], n(0); k<= kmaxwav[i]; ++k, ++n)
            {
                // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
                mu_wav_adapted_intinfo[i][n] = intinfo[i];
                switch (intinfo[i])
                {
                    case 0:
                        // LL: check whether mu is a left boundary wav
                        if (!(k < basis_->get_numofbw() ) ) 
                        {
                            // mui is not extended left
                            mu_wav_adapted_intinfo[i][n] = 4; // MM
                        }
                        break;
                    case 2:
                        // LR: check whether mui is a right boundary wav
                        if (!( (k > (1<<mu_j[i])- 1 - basis_->get_numofbw()) ) )
                        {
                            mu_wav_adapted_intinfo[i][n] = 1; // LM
                        }
                        break;
                    case 3:
                        // ML : check whether mu is a left boundary wav
                        if (!( k < basis_->get_numofbw()) )
                        {
                            // mui is not extended left. There is no intersection of mu and nu! This should not happen at this point
                            cout << "lambdanum = " << lambdanum 
                                    << "; lambda = " << *basis_->get_wavelet(lambdanum) 
                                    << "; levelnum = " << levelnum 
                                    << "; mu_j = " << mu_j << "; mu_p = " << mu_p << endl;
                            cout << "i = " << i << ";k = " << k << "; (1<<mu_j[i])- 1 - basis_->get_numofbw() = " << (1<<mu_j[i])- 1 - basis_->get_numofbw() << endl;
                            cout << "mu_wav_adapted_intinfo = " << mu_wav_adapted_intinfo 
                                    << "; mu_gen_adapted_intinfo = " << mu_gen_adapted_intinfo
                                    << "; intinfo = " << intinfo << endl;
                                    cout << "kminwav = " << kminwav << "; kmingen = " << kmingen << endl;
                            assert(false);
                        }
                        break;
                    case 5:
                        // MR: check whether mu is a right boundary wav
                        if (!(k > (1<<mu_j[i])- 1 - basis_->get_numofbw()) )
                        {
                            // nu is not extended right. There is no intersection of mu and nu! This should not happen at this point
                            
                            cout << "lambdanum = " << lambdanum 
                                    << "; lambda = " << *basis_->get_wavelet(lambdanum) 
                                    << "; levelnum = " << levelnum 
                                    << "; mu_j = " << mu_j << "; mu_p = " << mu_p << endl;
                            cout << "i = " << i << ";k = " << k << "; (1<<mu_j[i])- 1 - basis_->get_numofbw() = " << (1<<mu_j[i])- 1 - basis_->get_numofbw() << endl;
                            cout << "mu_wav_adapted_intinfo = " << mu_wav_adapted_intinfo 
                                    << "; mu_gen_adapted_intinfo = " << mu_gen_adapted_intinfo
                                    << "; intinfo = " << intinfo << endl;
                                    cout << "kminwav = " << kminwav << "; kmingen = " << kmingen << endl;
                            assert(false);
                        }
                        break;
                    case 6:
                        // RL : check whether mu is a left boundary gen/wav
                        if (!(k < basis_->get_numofbw()) )
                        {
                            // nu is not extended left. There is no intersection of mu and nu! This should not happen at this point
                            mu_wav_adapted_intinfo[i][n] = 7; // RM
                        }
                        break;
                    case 8:
                        // RR: check whether mu is a right boundary gen/wav
                        if (!(k > (1<<mu_j[i])- 1 - basis_->get_numofbw()) )
                        {
                            // nu is not extended right
                            mu_wav_adapted_intinfo[i][n] = 4; // MM
                        }
                        break;
                    case 1:
                    case 4:
                    case 7:
                        break;
                    default:
                        abort();
                        break;
                }
            } // end of loop over waveletBlock
        }
        // part 2b: use mu-adapted LMR-info and 1d integrals in wavelet/generatorBlocks to efficiently compute entries of w
        FixedArray1D< typename Block::iterator, DIM> block_it; //, block_it_begin;
        FixedArray1D< typename Array1D<unsigned int>::const_iterator ,DIM> mu_adapted_intinfo_it;
        typedef typename Index::type_type type_type;
        type_type current_type;
        for (unsigned int i=0;i<DIM;i++)
        {
            if (mu_min_type[i] == true)
            {
                // current_type[i] is coded as an integer, 0 = gen, 1 = wav
                current_type[i] = 0;
            }
            else
            {
                current_type[i] = 1;
            }
        }
        unsigned int number;
        Index temp_ind(basis_->first_wavelet(basis_->get_levelnum(mu_j,mu_p)));
        number = temp_ind.number(); // first wavelet with current_j,current_p
        FixedArray1D<int,DIM> jump_before, jump_after;
        unsigned int cjd; // current jump dimension == dimension (array index) where we are increasing k
        int blocksize;
        bool skip_this_type(false);
        for (unsigned int i=0; i<DIM; ++i)
        {
            if (mu_min_type[i] && (generatorBlock[i]->size() == 0))
            {
                skip_this_type=true;
                break;
            }
        }
        bool done = false;
        while (!done)
        {
            //compose the information in intersections into QTBasis::Index information and keep track of the number of the indices
            // for the current type e:
            // we iterate over all possible k
            // On the lowest dimension this is simple: every increase in k[DIM-1] increses the number of the wavelet by 1.
            // In other dimensions we need to keep track how the number changes, if we increase k[cjd] and set k[j] to the lowest possible value.
            // this is stored in jump_before and jump_after
            if (skip_this_type)
            {
                done = true; // == increase type e
                blocksize = 1;
                for (unsigned int i=0; i<DIM; ++i)
                {
                    if (current_type[i] == 0)
                    {
                        blocksize *= basis_->bases()[mu_p][i]->Deltasize(mu_j[i]);
                    }
                    else
                    {
                        blocksize *= basis_->bases()[mu_p][i]->Nablasize(mu_j[i]);
                    }
                }
                number += blocksize;
            }
            else
            {
                cjd=0;
                if (DIM == 1)
                {
                    if (current_type[0] == 0)
                    {
                        number += kmingen[0] - basis_->bases()[mu_p][0]->DeltaLmin();
                        blocksize = kmaxgen[0] - kmingen[0] +1;
                    } 
                    else
                    {
                        number += kminwav[0] - basis_->bases()[mu_p][0]->Nablamin();
                        blocksize = kmaxwav[0] - kminwav[0] +1;
                    }
                    for (unsigned int n = number; n < number + blocksize; ++n)
                    {
    // CLEANUP                        
                        assert (basis_->get_wavelet(n)->p() == mu_p);
                        assert (basis_->get_wavelet(n)->j() == mu_j);
                        assert (basis_->get_wavelet(n)->e() == current_type);
                        assert (basis_->get_wavelet(n)->k()[0] == ((current_type[0] == 0)?kmingen[0]:kminwav[0])+n-number);
                        assert (basis_->get_wavelet(n)->number() == n);
                        cout << "sorry! qtbasis::add_ball not completely implemented for 1d case" << endl;
                        abort();
                    }
                    assert ( number + ((current_type[0] == 0)? (basis_->bases()[mu_p][0]->DeltaRmax(mu_j[0]) - kmaxgen[0]) : (basis_->bases()[mu_p][0]->Nablamax(mu_j[0]) - kmaxwav[0]) )
                            ==
                            ((basis_->bases()[mu_p][0]->Deltasize(mu_j[0])) + ((current_type[0] == 1) ? (basis_->bases()[mu_p][0]->Nablasize(mu_j[0])):0)) );
                    number = ((basis_->bases()[mu_p][0]->Deltasize(mu_j[0])) + ((current_type[0] == 1) ? (basis_->bases()[mu_p][0]->Nablasize(mu_j[0])):0));
                    // if e[0] == 0: number now points to the first wavelet with the type e=1 but with the same (p,j). we expect that ++it does not increase p
                    // if e[0] == 1: this patch is finished and we expect that ++it increases p. Iteration should stop now
                    done = true; // i.e.,  increase type e at the end of the while loop
                }
                else
                {
                    int basf(0); // "blocks added so far" (on the current type e) via modulo calculus we can deduce which block we have to add next, i.e., where we have to increase k
                    while (!done)
                    {
                        if (cjd == 0)
                        {
                            // nothing has been done so far for this type
                            // compute all jump_numbers that are needed for the current type e
                            int temp_int(1);
                            for (int i = DIM-1; i >= 0; i--)
                            {
                                if (current_type[i] == 0)
                                { 
                                    jump_before[i] = temp_int;
                                    jump_before[i] *= kmingen[i] - basis_->bases()[mu_p][i]->DeltaLmin();
                                    jump_after[i] = temp_int;
                                    jump_after[i] *= basis_->bases()[mu_p][i]->DeltaRmax(mu_j[i]) - kmaxgen[i];
                                    temp_int *= basis_->bases()[mu_p][i]->Deltasize(mu_j[i]);
                                    block_it[i] = generatorBlock[i]->begin();
                                    mu_adapted_intinfo_it[i] = mu_gen_adapted_intinfo[i].begin();
                                }
                                else
                                {
                                    jump_before[i] = temp_int;
                                    jump_before[i] *= kminwav[i] - basis_->bases()[mu_p][i]->Nablamin();
                                    jump_after[i] = temp_int;
                                    jump_after[i] *= basis_->bases()[mu_p][i]->Nablamax(mu_j[i]) - kmaxwav[i];
                                    temp_int *= basis_->bases()[mu_p][i]->Nablasize(mu_j[i]);
                                    block_it[i] = waveletBlock[i]->begin();
                                    mu_adapted_intinfo_it[i] = mu_wav_adapted_intinfo[i].begin();
                                }
                            }
                            for (int i = DIM-2; i >= 0; --i)
                            {
                                jump_before[i] += jump_before[i+1];
                                jump_after[i] += jump_after[i+1];
                            }
                            if (current_type[DIM-1] == 0)
                            {
                                blocksize = kmaxgen[DIM-1] - kmingen[DIM-1] +1;
                            }
                            else
                            {
                                blocksize = kmaxwav[DIM-1] - kminwav[DIM-1] +1;
                            }
                            basf = 0; // no blocks with the current type were added so far
                        }
                        // "jump_before" for all dimensions (current_jump_dim, current_jump_dim+1,...,dim)
                        number += jump_before[cjd];
                        // add wavelets
                        // block_it points to the first_wavelet we want to add
    // CLEANUP                        
                        for (unsigned int n = number; n < number + blocksize; ++n)
                        {
                            assert (basis_->get_wavelet(n)->p() == mu_p);
                            assert (basis_->get_wavelet(n)->j() == mu_j);
                            assert (basis_->get_wavelet(n)->e() == current_type);
                            assert (basis_->get_wavelet(n)->k()[DIM-1] == ( (current_type[DIM-1] == 0)? kmingen[DIM-1]:kminwav[DIM-1]) +n-number); // this only checks the last entry of k
                            assert (basis_->get_wavelet(n)->number() == n);
                        }
    // TODO                        
                        // add wavelets
                        // compute_sum_over_patches (lambda_p, intinfo, integralshares);
                        compose_wavelets(w, number, blocksize, factor, lambda_p, 
                                mu_adapted_intinfo_it,
                                block_it,
                                precond);

                        
                        ++basf; // we have added a block!
                        number += blocksize;
                        // "increase k"
                        // find the first position (from the back) such that k[i] < kmax[i]
                        // then increase number by jump_after[j]. 
                        // Effect: the new number is the number of the first wavelet with the new value of k[i] and k[j]=kmin[j], j=i+1,..,DIM-1
                        // have we arrived at the very last possible k? -> done = true

                        // compute the dimension in which we want to increase k
                        int temp_i(basf); //, temp_count;
                        for (int i(DIM-1); i >= 0; --i)
                        {
                            if (i == 0) //we have already added the last possible block for this type, i.e., k[0] = maxk[0]
                            {
                                done = true; // apparently we have just added the last possible k for the current type vector e. We need to increase it
                                cjd = 0;
                                break;
                            }
                            if (current_type[i] == 0)
                            {
                                block_it[i] = generatorBlock[i]->begin();
                                mu_adapted_intinfo_it[i] = mu_gen_adapted_intinfo[i].begin();
                            }
                            else
                            {
                                block_it[i] = waveletBlock[i]->begin();
                                mu_adapted_intinfo_it[i] = mu_wav_adapted_intinfo[i].begin();
                            }
                            // use cjd as temporary int variable
                            cjd = (current_type[i-1] == 0)? (kmaxgen[i-1] - kmingen[i-1] +1):(kmaxwav[i-1]-kminwav[i-1]+1); // blocksize in this dimension
                            if ((temp_i % cjd) != 0)
                            {
                                cjd = i; // == "we can increase k[cjd]"
                                ++block_it[i-1];
                                ++mu_adapted_intinfo_it[i-1];
                                break;
                            }
                            temp_i = temp_i / cjd;
                        }
                        number += jump_after[cjd]; // number is the number of the first wavelet with the new value for k[cjd] and k[j] = Nablamin/Deltamin for j>cjd
                    }
                } // end of if (DIM > 1)
            } // end of   if (skip_this_type) {} else {}
            // increase the type e
            // number should be the number of the first wavelet with the new type
            // "small loop":
            // try to increase currenttype
            // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
            // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
            
            skip_this_type = false;
            int i(DIM-1);
            for (; i >= 0; i--)
            {
                // find first position on level j0
                if (mu_min_type[i])
                {
                    if (current_type[i] == 1)
                    {
                        current_type[i]=0;
                        if (generatorBlock[i]->size() == 0)
                        {
                            skip_this_type = true;
                        }
                        else
                        {
                            block_it[i] = generatorBlock[i]->begin();
                            mu_adapted_intinfo_it[i] = mu_gen_adapted_intinfo[i].begin();
                        }
                    }
// CLEANUP                    
                     /*
                    if (mu_min_type[DIM-1])
                    {
                        current_type[i] = 0;
                        block_it_begin[i] = generatorBlock[i]->begin();
                        mu_adapted_intinfo2[i] = & mu_gen_adapted_intinfo[i];
                    }
                    else
                    {
                        current_type[i] = 1;
                        block_it_begin[i] = waveletBlock[i]->begin();
                        mu_adapted_intinfo2[i] = & mu_wav_adapted_intinfo[i];
                    }
*/
                    else
                    {
                        //block_it_begin[i] = waveletBlock[i]->begin(); // value for i=0 is not used
                        current_type[i]=1;
                        block_it[i] = waveletBlock[i]->begin();
                        mu_adapted_intinfo_it[i] = mu_wav_adapted_intinfo[i].begin();
                        done = false;
                        break;
                    }
                }
            }
// CLEANUP
//            for (unsigned int k = 0; k < DIM; ++k)
//            {
//                cout << "mu_min_type[" << k << "] = " << mu_min_type[k] << "; current_type[" << k << "] = " << current_type[k] << "; generatorBlock[" << k << "]->size() = " << generatorBlock[k]->size() << endl;
//            }
            for (int j=0; j<i; ++j)
            {
                if ((current_type[j] == 0) && (generatorBlock[j]->size() == 0))
                {
                    skip_this_type=true;
                    break;
                }
            }
            // done == true means that all components with currentlevel[i]=j0[i] were wavelets, i.e.,
            // we have arrived at the type vector (1,1,...,1) and the iteration is finished
        } // end of while(true), i.e., all intersecting wavelets have been composed and added to the cache
    }
    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT, unsigned int DER_ONEDIMHAARCOUNT>
    void 
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT,DER_ONEDIMHAARCOUNT>
    ::compose_wavelets(Vector<double>& w,
                const unsigned int start,
                const unsigned int blocksize,
                const double factor,
                const unsigned int nu_p,
                const FixedArray1D< typename Array1D<unsigned int>::const_iterator ,DIM> mu_adapted_intinfo_it,
                const FixedArray1D< typename Block::iterator, DIM> block_it,
                const bool precond) const
                //const MultiIndex<unsigned int, DIM> intinfo, 
                //const FixedArray1D<entries,DIM> integralshares) const
    {
        assert ((DIM == 2) || (DIM == 3));
        assert (blocksize > 0);
        // The first OR the last wavelet in the block may have a different intinfo than the others.
        // Therefore, in this method, the geometry has to be computed exactly once or twice.
        unsigned int geometry_type;
        int centerpatchnumber;
        MultiIndex<bool, DIM> orientation;
        MultiIndex<unsigned int, DIM> intinfo;
        
        unsigned int current_offset(0);
        typename Block::iterator current_shares(block_it[DIM-1]);
        typename Array1D<unsigned int>::const_iterator current_last_intinfo (mu_adapted_intinfo_it[DIM-1]);
        bool recompute_geometry(false);
        unsigned int temp_i;
        double temp_d, temp_d1, temp_d2, temp_d3;
        while (current_offset < blocksize)
        {
            for (unsigned int i=0; i<DIM-1; ++i)
            {
                intinfo[i] = *mu_adapted_intinfo_it[i];
            }
            intinfo[DIM-1] = *current_last_intinfo;
            basis_->get_intersection_geometry(nu_p, intinfo, geometry_type, centerpatchnumber, orientation);

            if (DIM == 2)
            {
                int north, east, northeast;
                if (geometry_type == 3)
                {
                    if (centerpatchnumber == -1)
                    {
                        north = basis_->get_neighbours(nu_p,0);
                        east = basis_->get_neighbours(nu_p,2);
                        geometry_type = 8;
                    }
                    else
                    {
                        north = basis_->get_neighbours(centerpatchnumber,3);
                        east = basis_->get_neighbours(centerpatchnumber,1);
                        if (north != -1)
                        {
                            northeast = basis_->get_neighbours(north,1);
                            if (northeast == -1)
                            {
                                geometry_type = 10;
                            }
                            else
                            {
                                if (east == -1)
                                {
                                    geometry_type = 9;
                                }
                            }
                        }
                        else
                        {
                            northeast = basis_->get_neighbours(east,3);
                            geometry_type = 11;
                        }
                    }
                }
                switch (geometry_type)
                {
                    case 0: // 1 patch
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a_times_gram, rest; // <- this can be recycled for the most part of the block!
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a_times_gram[eta] = 0;
                            rest[eta] = 0;
                        }
                        if (orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                assert (a_times_gram[y] == 0);
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    a_times_gram[y] += agencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).first[x];
                                    rest[y] += qgencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).first[x]
                                            + agencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    a_times_gram[y] += agencoeffs_[centerpatchnumber](x,y) * (block_it[0]->second).first[ONEDIMHAARCOUNT-1-x];
                                    rest[y] += qgencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).first[ONEDIMHAARCOUNT-1-x]
                                            + agencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).second[ONEDIMHAARCOUNT-1-x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += a_times_gram[y] * (current_shares->second).second[y];
                                    temp_d += rest[y] * (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += a_times_gram[y] * (current_shares->second).second[ONEDIMHAARCOUNT-1-y];
                                    temp_d += rest[y] * (current_shares->second).first[ONEDIMHAARCOUNT-1-y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    case 1: // 2 patches: centerpatch and the patch right of it
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a_times_gram, rest;
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a_times_gram[eta] = 0;
                            rest[eta] = 0;
                        }
                        const int rightneighbour(basis_->get_neighbours(centerpatchnumber,1));
                        assert (rightneighbour != -1); // assert existence
                        if (!orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[rightneighbour] (x,y)) * (block_it[0]->second).first[x];
                                    rest[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[rightneighbour] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[rightneighbour] (x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[rightneighbour] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[rightneighbour] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[rightneighbour] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += a_times_gram[y] * (current_shares->second).second[y];
                                    temp_d += rest[y] * (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += a_times_gram[y] * (current_shares->second).second[ONEDIMHAARCOUNT-1-y];
                                    temp_d += rest[y] * (current_shares->second).first[ONEDIMHAARCOUNT-1-y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    case 2: // 2 patches: centerpatch and the one above
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a1_times_gram[eta] = 0;
                            a2_times_gram[eta] = 0;
                            rest1[eta] = 0;
                            rest2[eta] = 0;
                        }
                        const int north(basis_->get_neighbours(centerpatchnumber,3));
                        assert (north != -1); // assert existence
                        if (orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    a1_times_gram[y] += agencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).first[x];
                                    rest1[y] += qgencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).first[x]
                                             + agencoeffs_[centerpatchnumber] (x,y) * (block_it[0]->second).second[x];

                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a2_times_gram[y] += agencoeffs_[north] (x,y) * (block_it[0]->second).first[x];
                                    rest2[y] += qgencoeffs_[north] (x,y) * (block_it[0]->second).first[x]
                                             + agencoeffs_[north] (x,y) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    a1_times_gram[y] += agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                    rest1[y] += qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x]
                                             + agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).second[x];

                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a2_times_gram[y] += agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x];
                                    rest2[y] += qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).first[x]
                                             + agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[y] + rest2[ONEDIMHAARCOUNT-1-y])* (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[ONEDIMHAARCOUNT-1-y] + rest2[y])* (current_shares->second).first[y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    case 3: // 4 patches in a square. centerpatch is the lower left one
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a1_times_gram[eta] = 0;
                            a2_times_gram[eta] = 0;
                            rest1[eta] = 0;
                            rest2[eta] = 0;
                        }
                        if (orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (x,y) + qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[east] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[y] + rest2[ONEDIMHAARCOUNT-1-y])* (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[ONEDIMHAARCOUNT-1-y] + rest2[y])* (current_shares->second).first[y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    case 8: // 3 patches, L-shaped support. Southwest is missing, i.e., centerepatch == -1
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a1_times_gram[eta] = 0;
                            a2_times_gram[eta] = 0;
                            rest1[eta] = 0;
                            rest2[eta] = 0;
                        }
                        if (orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (x,y) + qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[east] (x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[east] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[east] (x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[y] + rest2[ONEDIMHAARCOUNT-1-y])* (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[ONEDIMHAARCOUNT-1-y] + rest2[y])* (current_shares->second).first[y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    case 9: // 3 patches, L-shaped support. Southeast is missing, i.e., east == -1
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a1_times_gram[eta] = 0;
                            a2_times_gram[eta] = 0;
                            rest1[eta] = 0;
                            rest2[eta] = 0;
                        }
                        if (orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (x,y) + qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (x,y) + agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[y] + rest2[ONEDIMHAARCOUNT-1-y])* (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[ONEDIMHAARCOUNT-1-y] + rest2[y])* (current_shares->second).first[y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    case 10: // 3 patches, L-shaped support. Northeast is missing, i.e., northeast == -1
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a1_times_gram[eta] = 0;
                            a2_times_gram[eta] = 0;
                            rest1[eta] = 0;
                            rest2[eta] = 0;
                        }
                        if (orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[east] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[north] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[northeast] (x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[y] + rest2[ONEDIMHAARCOUNT-1-y])* (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[ONEDIMHAARCOUNT-1-y] + rest2[y])* (current_shares->second).first[y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    case 11: // 3 patches, L-shaped support. North is missing, i.e., north == -1
                    {
                        FixedArray1D<double,ONEDIMHAARCOUNT> a1_times_gram, a2_times_gram, rest1, rest2;
                        for (unsigned int eta = 0; eta < ONEDIMHAARCOUNT; ++eta)
                        {
                            a1_times_gram[eta] = 0;
                            a2_times_gram[eta] = 0;
                            rest1[eta] = 0;
                            rest2[eta] = 0;
                        }
                        if (orientation[0])
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (x,y) + qgencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (x,y) + agencoeffs_[east] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[northeast] (ONEDIMHAARCOUNT-1-x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        else
                        {
                            for (int y = 0; y< ONEDIMHAARCOUNT; y++)
                            {
                                for (int x = 0 ; x<ONEDIMHAARCOUNT; x++)
                                {
                                    // agencoeffs_[centerpatchnumber] (i,j) == ith in x direction, jth in y direction!
                                    a1_times_gram[y] += (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * (block_it[0]->second).first[x];
                                    rest1[y] += (qgencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + qgencoeffs_[east] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[centerpatchnumber] (ONEDIMHAARCOUNT-1-x,y) + agencoeffs_[east] (x,y)) * (block_it[0]->second).second[x];
                                    a2_times_gram[y] += (agencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x];
                                    rest2[y] += (qgencoeffs_[northeast] (x,y)) * (block_it[0]->second).first[x]
                                            + (agencoeffs_[northeast] (x,y)) * (block_it[0]->second).second[x];
                                }
                            }
                        }
                        while(current_offset < blocksize)
                        {
                            temp_d = 0;
                            if (orientation[1])
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[y] + a2_times_gram[ONEDIMHAARCOUNT-1-y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[y] + rest2[ONEDIMHAARCOUNT-1-y])* (current_shares->second).first[y];
                                }
                            }
                            else
                            {
                                for (unsigned int y=0; y<ONEDIMHAARCOUNT; ++y)
                                {
                                    temp_d += (a1_times_gram[ONEDIMHAARCOUNT-1-y] + a2_times_gram[y]) * (current_shares->second).second[y];
                                    temp_d += (rest1[ONEDIMHAARCOUNT-1-y] + rest2[y])* (current_shares->second).first[y];
                                }
                            }
                            w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                            temp_i = *current_last_intinfo;
                            ++current_last_intinfo;
                            ++current_shares;
                            ++current_offset;
                            if (current_offset != blocksize)
                            {
                                if (temp_i != *current_last_intinfo)
                                {
                                    assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                                    recompute_geometry = true;
                                    break;
                                }
                            }
                        }
                    }
                        break;
                    default:
                        cout << "this geometry is not implemented" << endl;
                        abort();
                        break;
                } // end of switch(geometry))
            } // end of DIM == 2
            else
            {
                // DIM == 3. Only Poisson equation implemented!
                if (ONEDIMHAARCOUNT != 1)
                {
                    cout << "warning: calling compute_wavelets with ONEDIMHAARCOUNT > 1 is inefficient. \nThis method computes the Poisson equation, so ==1 is sufficient!" << endl;
                    abort();
                }
                switch (geometry_type)
                {
                    case 0: // 1 patch
                        temp_d3 = 1;
                        break;
                    case 1:
                    case 2:
                    case 4: // 2 patches
                        temp_d3 = 2;
                        break;
                    case 3:
                    case 5:
                    case 6: // 4 patches
                        temp_d3 = 4;
                        break;
                    case 7: // 8 patches
                        temp_d3 = 8;
                        break;
                    default:
                        abort();
                        break;
                }
                
                temp_d1 = (block_it[0]->second).second[0] * (block_it[1]->second).first[0]
                        + (block_it[0]->second).first[0] * (block_it[1]->second).second[0];
                temp_d2 = (block_it[0]->second).first[0] * (block_it[1]->second).first[0];
                while(current_offset < blocksize)
                {
                    temp_d = temp_d1 * (current_shares->second).first[0]
                            + temp_d2 * (current_shares->second).second[0];
                    temp_d *= temp_d3;

                    w[start+current_offset] += precond? (factor*temp_d/D(start+current_offset)) : (factor*temp_d);
                    temp_i = *current_last_intinfo;
                    ++current_last_intinfo;
                    ++current_shares;
                    ++current_offset;
                    if (current_offset != blocksize)
                    {
                        if (temp_i != *current_last_intinfo)
                        {
                            assert (recompute_geometry == false); // otherwise we recompute the geometry for the 2nd time - this should not happen since Block[DIM-1] contains only boundary wavelets at the start OR at the end!
                            recompute_geometry = true;
                            break;
                        }
                    }
                }
            } // end of else DIM == 3; geometry will be recomputed if current_offset < blocksize
        } // end of while (current_offset < blocksize)
    }
    
    
// ------------------------------------------------------------


#if 0


    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT>
    void 
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT>::evaluate(Array1D<SampledMapping<QTBASIS::space_dimension> > & previous_sampling,
            const unsigned int lambdanum,
            const bool primal,
            const int resolution)
    {
        
    }
        
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT>
    Array1D<SampledMapping<QTBASIS::space_dimension> > 
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT>::evaluate(const InfiniteVector<double, int>& coeffs,
            const bool primal,
            const int resolution)
    {
        /*
        typedef typename TensorBasis<IBASIS,DIM>::Index Index;
        for (typename InfiniteVector<double,Index>::const_iterator it(coeffs.begin()),
               itend(coeffs.end()); it != itend; ++it)
               {
                   result.add(*it, evaluate(basis, it.index(), primal, resolution));
               }
         */
    }

#endif
     

#if 0
        // old code from a(,)
        // computes entries of a (!!)DIM(!!)-dimensional cache as a tensor product of given one dimensional index sets.
        typedef std::map<int, double> Block;
        typedef std::map<int, Block> Column;
        typedef std::map<int, Column> ColumnCache;
        ColumnCache entries_cache_;
        
        typename ColumnCache::iterator col_lb(entries_cache_.lower_bound(nunum));
        typename ColumnCache::iterator col_it(col_lb);

        if (col_lb == entries_cache_.end() ||
            entries_cache_.key_comp()(nunum, col_lb->first))
        {
            // this column has not been requested so far
            // insert a new column and continue with the blocks
            typedef typename ColumnCache::value_type value_type;
            col_it = entries_cache_.insert(col_lb, value_type(nunum, Column()));
        }
        
        
        // Be careful, there is no generator level in the cache! The situation from the MRA setting:
        // KEY OF GENERATOR LEVEL IS j0-1 NOT j0 !!!!
        // does not hold in the tensor setting. Generators and wavelets on
        // the minimal level are thrown together in one index set (componentwise),
        // this also applies to tensors of generators on the lowest level
        
        //const int lambda_num = lambda.number();
        //const int nu_num = nu.number();


        // search for column 'mu'
        typename ColumnCache::iterator col_lb(entries_cache_.lower_bound(nunum));
        typename ColumnCache::iterator col_it(col_lb);

        if (col_lb == entries_cache_.end() ||
            entries_cache_.key_comp()(nunum, col_lb->first))
        {
            // this column has not been requested so far
            // insert a new column and continue with the blocks
            typedef typename ColumnCache::value_type value_type;
            col_it = entries_cache_.insert(col_lb, value_type(nunum, Column()));
        }

// CLEANUP        
/*        
        index_lt lambda_key(lambda.j()),first_level(basis().j0());
        for (int k=0;k<space_dimension;k++)
        {
            lambda_key[k] = lambda_key[k]-first_level[k];
        }
//TODO (PERFORMANCE): store the numbers of all levels up to jmax, do not compute anything here:
        int blocknumber(lambda_key.number());
 * 
 * 
        // check wether entry has already been computed
        typedef std::list<Index> IntersectingList;
*/
        
        

        // col_it points to the column of \mu
        Column& col(col_it->second);
        // check whether the level 'lambda' belongs to has already been calculated
        Index nu(basis_->get_wavelet(nunum));
        int levelblocknumber(basis_->get_levelnum(nu.j(),nu.p()));
        typename Column::iterator lb(col.lower_bound(levelblocknumber));
        typename Column::iterator it(lb);

        if (lb == col.end() ||
            col.key_comp()(levelblocknumber, lb->first))
        {
            // no entries have ever been computed for this column and this level
            // insert a new block, 
            // then compute whole level block (a(mu,nu))_mu for all mu that intersect nu. 
            // The entries of A are decomposable into sums over Haar wavelet coefficients, since A is determinded by acoeffs_ and qcoeffs_.
            // For each fixed Haar wavelet the integrals over all mu in the Block are decomposable into one dimensional integrals.
            // These one dimensional integrals have their own cache, since they are reused several times.
            // The code has some similarities with intersecting_wavelets from tbasis_support.
            typedef typename Column::value_type value_type;
            it = col.insert(lb, value_type(levelblocknumber, Block()));
            Block& block(it->second);

            
            
// CLEANUP
            //for the relevant patches: for the relevant types e_ : a pair (k_first,k_last) consisting of the first and last k that intersects
            //the intersecting wavelets are all wavelets with k1 <= k <= k2
            //typedef typename std::list<std:: pair <unsigned int, std::pair< MultiIndex<int,DIM>, std::pair< MultiIndex<int,DIM>, MultiIndex<int,DIM> > > > > IntersectionList;
            //list<std:: pair <unsigned int, std::pair< MultiIndex<int,DIM>, std::pair< MultiIndex<int,DIM>, MultiIndex<int,DIM> > > > > 
            //IntersectionList intersections; 
          
            
            
            // For each fixed type e_, the set of indices mu such that supp psi_mu intersects supp psi_nu 
            // is a cartesian product of sets corresponding to one dimensional support intersections.
            // Depending on e[i] these one dimensional intersecting sets are fully characterized by the first and last gen/wav in it.
            MultiIndex<int,DIM> kmingen, kminwav, kmaxgen, kmaxwav;
            MultiIndex<unsigned int, DIM> intinfo;
// CLEANUP            
//            basis_->intersecting_wavelets(nunum, levelblocknumber, kmingen, kmaxgen, kminwav, kmaxwav, intinfo);
//            
// CHECK USAGE OF lambda_p, mu_p, min_type [which patch?]            
            
            // insert all intersecting wavelets into the cache (compare with tbasis_support::intersecting_wavelets)
            // This is an iteration over all possible types e_ of the mu's
            
// CLEANUP

            //typename IntersectionList::const_iterator it(intersections.begin()), itend(intersections.end());
            //assert (it != itend);
            
            
            
            //MultiIndex<int,DIM> j(lambda.j()), e((*it).first), kmin((*it).second.first), kmax((*it).second.second);
            //unsigned int p(*it); 
            
            
            MultiIndex<int,DIM> j(nu.j());
            unsigned int p(nu.p());
            bool generators = true;
            typedef typename Index::type_type type_type;
            type_type current_type;
            MultiIndex<bool,DIM> min_type;
            for (unsigned int i=0;i<DIM;i++)
            {
                min_type[i]=(j[i]==basis_->j0()[p][i]); // remember whether we are on the minimal levels to avoid calls of basis_->j0()
                current_type[i]=(min_type[i] == true) ? 0:1; // current_type[i] is coded as an integer, 0 = gen, 1 = wav
                generators = generators && min_type[i];
            }
            unsigned int number;
            Index temp_ind(basis_->first_wavelet(basis_->get_levelnum(j,p)));
            number = temp_ind.number(); // first wavelet with current_j,current_p
            FixedArray1D<int,DIM> jump_before, jump_after;
            unsigned int cjd;
            int blocksize;

            while (true)
            {
                //compose the information in intersections into QTBasis::Index information and keep track of the number of the indices
                // for the current type e:
                // we iterate over all possible k
                // On the lowest dimension this is simple: every increase in k[DIM-1] increses the number of the wavelet by 1.
                // In other dimensions we need to keep track how the number changes, if we increase k[cjd] and set k[j] to the lowest possible value.
                // this is stored in jump_before and jump_after
                cjd=0;
                bool done = false;
                if (DIM == 1)
                {
                    if (current_type[0] == 0)
                    {
                        number += kmingen[0] - basis_->bases()[p][0]->DeltaLmin();
                        blocksize = kmaxgen[0] - kmingen[0] +1;
                    } else
                    {
                        number += kminwav[0] - basis_->bases()[p][0]->Nablamin();
                        blocksize = kmaxwav[0] - kminwav[0] +1;
                    }
                    for (unsigned int n = number; n < number + blocksize; ++n)
                    {
// CLEANUP                        
                        assert (basis_->get_wavelet(n).p() == p);
                        assert (basis_->get_wavelet(n).j() == j);
                        assert (basis_->get_wavelet(n).e() == current_type);
                        assert (basis_->get_wavelet(n).k() == ((current_type[0] == 0)?kmingen[0]:kminwav[0])+n-number);
                        assert (basis_->get_wavelet(n).number() == n);
                        abort();
// TODO                        
                        // add wavelets
                        //intersecting.push_back(basis.get_wavelet(n));
                    }
                    assert ( number + ((current_type[0] == 0)? (basis_->bases()[p][0]->DeltaRmax(j[0]) - kmaxgen[0]) : (basis_->bases()[p][0]->Nablamax(j[0]) - kmaxwav[0]) ) 
                            ==
                            ((basis_->bases()[p][0]->DeltaSize(j[0])) + ((current_type[0] == 1) ? (basis_->bases()[p][0]->NablaSize(j[0])):0)) );
                    number = ((basis_->bases()[p][0]->DeltaSize(j[0])) + ((current_type[0] == 1) ? (basis_->bases()[p][0]->NablaSize(j[0])):0));
                    // if e[0] == 0: number now points to the first wavelet with the type e=1 but with the same (p,j). we expect that ++it does not increase p
                    // if e[0] == 1: this patch is finished and we expect that ++it increases p. Iteration should stop now
                    
                    
                    done = true; // i.e.,  increase type e at the end of the while loop
                }
                else
                {
                    int basf(0); // "blocks added so far" (on the current type e) via modulo calculus we can deduce which block we have to add next, i.e., where we have to increase k
                    while (!done)
                    {
                        if (cjd == 0)
                        {
                            // nothing has been done so far for this type
                            // compute all jump_numbers that are needed for the current type e
                            int temp_int(1);
                            for (int i = DIM-1; i >= 0; i--)
                            {
                                if (current_type[i] == 0)
                                {
                                    jump_before[i] = temp_int;
                                    jump_before[i] *= kmingen[i] - basis_->bases()[p][i]->DeltaLmin();
                                    jump_after[i] = temp_int;
                                    jump_after[i] *= basis_->bases()[p][i]->DeltaRmax(j[i]) - kmaxgen[i];
                                    temp_int *= basis_->bases()[p][i]->Deltasize(j[i]);
                                }
                                else
                                {
                                    jump_before[i] = temp_int;
                                    jump_before[i] *= kminwav[i] - basis_->bases()[p][i]->Nablamin();
                                    jump_after[i] = temp_int;
                                    jump_after[i] *= basis_->bases()[p][i]->Nablamax(j[i]) - kmaxwav[i];
                                    temp_int *= basis_->bases()[p][i]->Nablasize(j[i]);
                                }
                            }
                            for (int i = DIM-2; i >= 0; --i)
                            {
                                jump_before[i] += jump_before[i+1];
                                jump_after[i] += jump_after[i+1];
                            }
                            if (current_type[DIM-1] == 0)
                            {
                                blocksize = kmaxgen[DIM-1] - kmingen[DIM-1] +1;
                            }
                            else
                            {
                                blocksize = kmaxwav[DIM-1] - kminwav[DIM-1] +1;
                            }
                            basf = 0; // no blocks with the current type were added so far
                        }
                        // "jump_before" for all dimensions (current_jump_dim, current_jump_dim+1,...,dim)
                        number += jump_before[cjd];
                        // add wavelets
                        for (unsigned int n = number; n < number + blocksize; ++n)
                        {
    // CLEANUP                        
                            assert (basis_->get_wavelet(n).p() == p);
                            assert (basis_->get_wavelet(n).j() == j);
                            assert (basis_->get_wavelet(n).e() == current_type);
                            assert (basis_->get_wavelet(n).k()[DIM-1] == ( (current_type[DIM-1] == 0)? kmingen[DIM-1]:kminwav[DIM-1]) +n-number); // this only checks the last entry of k
                            assert (basis_->get_wavelet(n).number() == n);
                            abort();
    // TODO                        
                            // add wavelets
                            //intersecting.push_back(basis.get_wavelet(n));
                            
                        }
                        ++basf; // we have added a block!
                        number += blocksize;
                        // "increase k"
                        // find the first position (from the back) such that k[i] < kmax[i]
                        // then increase number by jump_after[j]. 
                        // Effect: the new number is the number of the first wavelet with the new value of k[i] and k[j]=kmin[j], j=i+1,..,DIM-1
                        // have we arrived at the very last possible k? -> done = true

                        // compute the dimension in which we want to increase k
                        int temp_i(basf);
                        for (int i(DIM-1); i >= 0; --i)
                        {
                            if (i == 0) //we have already added the last possible block for this type, i.e., k[0] = maxk[0]
                            {
                                done = true; // apparently we have just added the last possible k for the current type vector e. We need to increase it
                                cjd = 0;
                                break;
                            }
                            int temp_count = (current_type[i-1] == 0)? (kmaxgen[i-1] - kmingen[i-1] +1):(kmaxwav[i-1]-kminwav[i-1]+1); // blocksize in this dimension
                            if ((temp_i % temp_count) != 0)
                            {
                                cjd = i; // == "we can increase k[cjd]"
                                break;
                            }
                            temp_i = temp_i / temp_count;
                        }
                        number += jump_after[cjd]; // number is the number of the first wavelet with the new value for k[cjd] and k[j] = Nablamin/Deltamin for j>cjd
                    }
                } // end of if (DIM > 1)
                // increase the type e
                // number should be the number of the first wavelet with the new type
                // "small loop":
                // try to increase currenttype
                // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
                // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
                for (int i(DIM-1); i >= 0; i--)
                {
                    // find first position on level j0
                    if (min_type[i])
                    {
                        if (current_type[i] == 1)
                        {
                            current_type[i]=0;
                        } else
                        {
                            current_type[i]=1;
                            done = false;
                            break;
                        }
                    }
                }
                // done == true means that all components with currentlevel[i]=j0[i] were wavelets, i.e.,
                // we have arrived at the type vector (1,1,...,1) and the iteration is finished
            } // end of while(true), i.e., all intersecting wavelets have been composed and added to the cache
/*            
            if ((blocknumber == 0) && (lambda.j() == first_level))
            {
                // Insert Generators and Wavelets

                // also add Generators to the Cache
                // this case has to be considered seperatly because
                // intersecting_wavelets divides the case of a basis
                // function made entirely of generators and the cache does not.

                // TODO : produce nicer looking code:
                IntersectingList nusG,nusW;
                intersecting_wavelets(basis(), nu,
                                      lambda.j(),
                                      true, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tbasis_support!
                                      nusG);
                intersecting_wavelets(basis(), nu,
                                      lambda.j(),
                                      false, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tbasis_support!
                                      nusW);
                // compute entries
                for (typename IntersectingList::const_iterator it(nusG.begin()), itend(nusG.end());it != itend; ++it)
                {
                    const double entry = problem->a(*it, nu);
                    typedef typename Block::value_type value_type_block;
                    if (fabs(entry) > 1e-16 ) //(entry != 0.)
                    {
                        block.insert(block.end(), value_type_block((*it).number(), entry));
                        if ((*it).number() == lambda_num)
                        {
                            r = entry;
                        }
                    }
                }

                for (typename IntersectingList::const_iterator it(nusW.begin()), itend(nusW.end());it != itend; ++it)
                {
                    const double entry = problem->a(*it, nu);
                    typedef typename Block::value_type value_type_block;
                    if (fabs(entry) > 1e-16 ) //(entry != 0.)
                    {
                        // Insertion should be efficient, since the wavelets coma after the generators
                        block.insert(block.end(), value_type_block((*it).number(), entry));
                        if ((*it).number() == lambda_num)
                        {
                            r = entry;
                        }
                    }
                }
            }
            else
            {
                // there are no Generators

                IntersectingList nus;
                intersecting_wavelets(basis(),
                                      nu,
                                      //std::max(j, basis().j0()),
                                      lambda.j(),
                                      false, // this argument isn't that helpful for the way the method is used in a() .  It doesn't play a role for DIM>1 and miserably fails for DIM=1. Need to work on the routine in tbasis_support!
                //                    == (basis().j0()-1),
                                      nus);
                // compute entries
                for (typename IntersectingList::const_iterator it(nus.begin()), itend(nus.end());it != itend; ++it)
                {
                    const double entry = problem->a(*it, nu);
                    typedef typename Block::value_type value_type_block;
                    if (fabs(entry) > 1e-16 ) //(entry != 0.)
                    {
                        block.insert(block.end(), value_type_block((*it).number(), entry));
                        if ((*it).number() == lambda_num)
                        {
                            r = entry;
                        }
                    }
                }
            }
        
*/        
        } // end of "insert new levelblock"
        
        
        
        // level already exists --> extract row corresponding to 'lambda'
        else
        {
            Block& block(it->second);

            //typename Block::iterator block_lb(block.lower_bound(lambda));
            typename Block::iterator block_lb(block.lower_bound(nunum));
            typename Block::iterator block_it(block_lb);
            // level exists, but in row 'lambda' no entry is available ==> entry must be zero
            if (block_lb == block.end() ||
                block.key_comp()(nunum, block_lb->first))
            {
                r = 0;
            }
            else 
            {
                r = block_it->second;
            }
        }
#endif
        
#if 0    
    template <class QTBASIS, unsigned int ONEDIMHAARCOUNT>
    void
    CachedQTProblem<QTBASIS,ONEDIMHAARCOUNT>::add_ball(const Index& lambda,
                                                   Vector<double>& w,
                                                   const int radius,
                                                   const double factor,
                                                   const int maxlevel,
                                                   const CompressionStrategy strategy,
                                                   const bool precond)
    {
        /*
         * For dim=1,2 it is easy to describe all levels in the (1-norm) ball of radius around lambda.
         * (Attention: possibly different minimal levels on different patches!)
         * For higher dimensions this is done recursivly.
         */
        int lambda_num = lambda.number();
        double d1 = precond ? D(lambda_num) // stored matrix is unpreconditioned
                            : 1.0;

        MultiIndex<int,DIM> temp_j;
        unsigned int temp_p;
        if (DIM == 1)
        {
            basis_->get_level(0, temp_j,temp_p);
            //unsigned int start_levelnum(0);
            unsigned int start_num(0);
            if (temp_j[0] < (lambda.j()[0]-radius))
            {
                temp_j[0] = (lambda.j()[0]-radius);
                //start_levelnum = basis_->get_levelnum(temp_j,0); // uses: minimal levels differ at most by 1
                start_num = basis_->first_wavelet(basis_->get_levelnum(temp_j,0)); // uses: minimal levels differ at most by 1
                //int start_index = (temp_j[0] >= (lambda.j()[0]-radius))? 0 : basis().first_wavelet(lambda.j()[0]-radius).number();
            }
            temp_j[0] = min(lambda.j()[0]+radius,basis_->get_jmax());
            //unsigned int end_levelnum(basis_->get_levelnum(temp_j, basis_->get_nop()-1));
            unsigned int end_num(basis_->last_wavelet( basis_->get_levelnum(temp_j, basis_->get_nop()-1)) );
            //int end_index = basis().last_wavelet(min(lambda.j()[0]+radius,basis_->get_jmax())).number();
            double entry;
            for (unsigned int i = start_num; i <= end_num; i++)
            {                
                entry = a(i,lambda_num); // this is too expensive!  But 1d case isn't important
                if (entry != 0)
                {
                    w[i] += precond ? (entry * factor / (d1 * D(i,i) ))
                                    : entry * factor;
                }
            }
        }
        else if (DIM == 2)
        {
            // The ball can be described of levellines consisting of levels with the same multidegree
            // The first is determined with the distance of lambda.j and (minimal) j0. The last with basis_->get_jmax()
            // The first level in a levelline is determined with minx = max(j0[0], lambda.j[0]-radius)
            // The last level in a levelline is determined with miny = max(j0[1], lambda.j[1]-radius)
            
            basis_->get_level(0, temp_j,temp_p);
            
            //index_lt j0(this->basis().j0());
            int lambdaline = lambda.j()[0]+lambda.j()[1];
            int lowestline = temp_j[0]+temp_j[1];
            int dist2j0=lambdaline-lowestline;
            int dist2maxlevel=basis_->get_jmax()-lambdaline;
            
            MultiIndex<int,DIM> min_j, max_j;
            int xstart,xend,ystart,steps;
            // iterate the levellines. offset relative to lambdas levelline
            //int start_num, end_num;
            unsigned int first_levelnum, last_levelnum; // first and last level on the current levelline
            for (int offset = -std::min(dist2j0,radius); offset < std::min(dist2maxlevel,radius)+1; offset++)
            {
                // iterate over the levels on the levelline

                // ignoring restrictions by j0 for the moment, we have:
                xstart = lambda.j()[0]-radius+ceil((radius+offset)/2.0); //x coordinate of the first level
                xend = lambda.j()[0]+floor((radius+offset)/2.0); // same for the last
                ystart = lambda.j()[1]+floor((radius+offset)/2.0); // and for the second dimension

                // we walk (make steps) through the levelline
                // The first level in the current levelline hast steps=0.
                // restrictions by j0 mean:
                // The first level on the current levelline may have steps >0.
// CLEANUP
                assert (temp_j[0] == basis_->j0()[patch_with_minimal_j0_[0].front()][0]);
                
                bool is_minimal_x(false), is_minimal_y(false);
                steps = max(0,temp_j[0]-xstart);
                min_j[0]= xstart+steps;
                min_j[1]= ystart-steps;
                if (temp_j[0] <= xstart)
                {
                    is_minimal_x = true;
                }
                steps = min(xend-xstart,ystart-basis_->j0()[patch_with_minimal_j0_[1].front()][1]);
                max_j[0]= xstart+steps;
                max_j[1]= ystart-steps;
                
                if (ystart-basis_->j0()[patch_with_minimal_j0_[1].front()][1] >= xend -xstart )
                {
                    is_minimal_y = true;
                }
                /*
                 * min_j, max_j contain the minimal and maximal level j on the current levelline
                 * iterate over all those levels. The first level may have the minimal j0[x-direction]. 
                 * In this case we iterate only over those patches with this minimal value.
                 * In the last step the analogous restriction in the y direction applies
                 */
                if (is_minimal_x)
                {
                    if (is_minimal_y)
                    {
                        // iterate over patches in the intersection of patch_with_minimal_j0_[0] and patch_with_minimal_j0_[1]
                        list<int> intersection;
                        set_intersection(patch_with_minimal_j0_[0].begin(), patch_with_minimal_j0_[0].end(),
                                patch_with_minimal_j0_[1].begin(), patch_with_minimal_j0_[1].end(),
                                intersection);
    // cleanup
                        for (list<int>::const_iterator it(intersection.begin()), it2(intersection.begin()), itend(intersection.end()); it!= itend; ++it)
                        {
                            ++it2;
                            if (it2 != itend)
                            {
                                assert (*it < *it2);
                            }
                        }
                        assert (min_j == max_j);
                        
                        first_levelnum = basis_->get_levelnum(min_j, *intersection.begin());
                        last_levelnum = basis_->get_levelnum(min_j, *intersection.end());
                    }
                    else
                    {
                        // iterate over patches in patch_with_minimal_j0_[0]
                        first_levelnum = basis_->get_levelnum(min_j, patch_with_minimal_j0_[0].begin());
                        last_levelnum = basis_->get_levelnum(max_j, basis_->get_nop()-1);
                    }
                }
                else
                {
                    if (is_minimal_y)
                    {
                        // iterate over the patches in minimal_level_y
                        first_levelnum = basis_->get_levelnum(min_j, basis_->get_nop()-1);
                        last_levelnum = basis_->get_levelnum(max_j, patch_with_minimal_j0_[1].end());
                    }
                    else
                    {
                        // iterate over all patches
                        first_levelnum = basis_->get_levelnum(min_j, basis_->get_nop()-1);
                        last_levelnum = basis_->get_levelnum(max_j, basis_->get_nop()-1);
                    }
                }
                
                for (unsigned int levelnum = first_levelnum; levelnum <= last_levelnum; levelnum++)
                {
                    add_level(lambda_num,
                            w,
                            levelnum,
                            factor/d1);
                }
// TODO: multiply w[i] with D[i]^{-1/2} 
                abort();
            }
        }
        else // dim > 2. iteration over all levels in 'range' is done recursive for arbitrary dimensions
        {
            //add_level_recurse(lambda,w,radius,factor,lambda.j(),0,basis_->get_jmax(),true,d1,precond);
            cout << "Cached_QTPROBLEM::add_ball is not jet implemented" << endl;
            cout << "see CachedTProblem::add_ball for details" << endl;
            abort();
        }
    }
#endif

    
}

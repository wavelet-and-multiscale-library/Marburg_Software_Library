
#include "qtbasis_index.h"

// implementation for qtbasis.h

namespace WaveletTL
{
    template <class IBASIS, unsigned int DIM>
    QTBasis<IBASIS,DIM>::QTBasis() : numofbw_((IBASIS::primal_polynomial_degree() + IBASIS::primal_vanishing_moments() -2) / 2)
    {
        // we only need one instance of IBASIS (with homogeneous b.c.)
    	IBASIS* b = new IBASIS(true,true);
    	bases_infact_[0] = b;
        //boundary_gens_infact_[0]=boundary_gens_infact_[1]=boundary_wavs_infact_[0]=boundary_wavs_infact_[1]=0; // no need to initialize _infact for the other bases
        //j0 = new Array1D<MultiIndex<int,DIM>> (1);
        j0_.resize(1);
        bases_.resize(1);
        corners_.resize(1); // origin should be set correctly by default Point<DIM> constructor
        neighbours_.resize(1);
        extinfo_.resize(1);
        //boundary_gens_.resize(1);
        //boundary_wavs_.resize(1);
        //j0_[0] = b->j0();
        //j0_[0] = new MultiIndex<int,DIM> ();
        //bases_[0] = new FixedArray1d<IBASIS*,DIM>
    	for (unsigned int i = 0; i < DIM; i++)
        {
            bases_[0][i] = b;
            j0_[0][i] = b->j0();
            neighbours_[0][2*i] = -1;
            neighbours_[0][2*i+1] = -1;
            extinfo_[0][2*i] = 0;
            extinfo_[0][2*i+1] = 0;
        }
        intersection_geometry_centerpatchnumber_.resize(1,intpower(9,DIM));
        intersection_geometry_type_.resize(1,intpower(9,DIM));
        for (unsigned int i=0; i<DIM; ++i)
        {
            intersection_geometry_orientation_[i].resize(1,intpower(9,DIM));
        }
        for (unsigned int j = 0; j<intpower(9,DIM); j++)
        {
            intersection_geometry_centerpatchnumber_.set_entry(0,j,-1);
        }
    }

    template <class IBASIS, unsigned int DIM>
    QTBasis<IBASIS,DIM>::QTBasis(const Array1D<Point<DIM, int> >& corners,
            const Array1D<FixedArray1D<int,2*DIM> >& neighbours,
            const Array1D<FixedArray1D<bool,2*DIM> >& bc)
    : corners_(corners), 
            neighbours_(neighbours), 
            numofbw_((IBASIS::primal_polynomial_degree() + IBASIS::primal_vanishing_moments() -2) / 2),
            bc_(bc)
    {
        assert (corners_.size() == neighbours_.size() && neighbours_.size() == bc.size());

        j0_.resize(corners_.size());
        bases_.resize(corners_.size());
        extinfo_.resize(corners_.size());

        assert (bc.size() == bc_.size());
        for (unsigned int i=0; i< bc_.size(); i++)
        {
            for (unsigned int j=0;j<2*DIM; ++j)
            {
                assert(bc[i][j] == bc_[i][j]);
            }
        }
        int tempint;
        IBASIS* b;
        for (unsigned int i=0; i< 4;++i)
        {
            bases_infact_[i]=0;
        }
        for (unsigned int p=0; p < corners_.size();++p)
        {
            for (unsigned int i = 0; i < DIM; i++)
            {
                tempint = ((bc[p][2*i])?0:2) + ((bc[p][2*i+1])?0:1);
                // check whether the corresponding 1d basis already exists
                b = bases_infact_[tempint];
                if (b==0)
                {
                    b = new IBASIS(bc[p][2*i], bc[p][2*i+1]);
                    bases_infact_[tempint] = b;
#if 0                    
                    // compute number of nonvanishing boundary generators/wavelets
                    if (bc[p][2*i] == true)
                    {
                        boundary_wavs_infact_[2*tempint] = boundary_gens_infact_[2*tempint] = boundary_wavs_[p][2*i] = boundary_gens_[p][2*i] = 0;
                    }
                    else
                    {
                        tempind = b->first_generator(b->j0());
                        while(true)
                        {
                            if (abs(b->evaluate(0,tempind,0.0)) > 1e-16)
                            {
                                ++boundary_gens_infact_[2*tempint];
                                ++tempind;
                            }
                            else
                            {
                                break;
                            }
                        }
                        assert (boundary_gens_infact_[2*tempint] > 0);
                        tempind = b->first_wavelet(b->j0());
                        while(true)
                        {
                            if (abs(b->evaluate(0,tempind,0.0)) > 1e-16)
                            {
                                ++boundary_wavs_infact_[2*tempint];
                                ++tempind;
                            }
                            else
                            {
                                break;
                            }
                        }
                        assert (boundary_wavs_infact_[2*tempint] > 0);
                    }
                    if (bc[p][2*i+1] == true)
                    {
                        boundary_wavs_infact_[2*tempint+1] = boundary_gens_infact_[2*tempint+1] = boundary_wavs_[p][2*i+1] = boundary_gens_[p][2*i+1] = 0;
                    }
                    else
                    {
                        tempind = b->first_generator(b->j0());
                        tempind2 = b->last_generator(b->j0());
// CLEANUP
                        //cout << (abs(b->evaluate(0,tempind,1.0))) << endl;
                        //cout << (abs(b->evaluate(0,tempind,1.0)) < 1e-16) << endl;
                        //cout << tempind << " != " << tempind2 << " = " << (tempind != tempind2) << endl;
                        while((abs(b->evaluate(0,tempind,1.0)) < 1e-16) && (tempind != tempind2))
                        {
                            ++tempind;
                        }
                        while(true)
                        {
                            ++boundary_gens_infact_[2*tempint+1];
                            if (tempind == tempind2)
                            {
                                break;
                            }
                            else
                            {
                                ++tempind;
                            }
                        }
                        assert (boundary_gens_infact_[2*tempint+1] > 0);
                        
                        tempind = b->first_wavelet(b->j0());
                        tempind2 = b->last_wavelet(b->j0());
                        while((abs(b->evaluate(0,tempind,1.0)) < 1e-16) && tempind != tempind2)
                        {
                            ++tempind;
                        }
                        while(true)
                        {
                            ++boundary_wavs_infact_[2*tempint+1];
                            if (tempind == tempind2)
                            {
                                break;
                            }
                            else
                            {
                                ++tempind;
                            }
                        }
                        assert (boundary_wavs_infact_[2*tempint+1] > 0);
                    }
#endif
                }
                bases_[p][i] = b;
                j0_[p][i] = b->j0();
                //boundary_gens_[p][2*i]=boundary_gens_infact_[2*tempint];
                //boundary_wavs_[p][2*i]=boundary_wavs_infact_[2*tempint];
                //boundary_gens_[p][2*i+1]=boundary_gens_infact_[2*tempint+1];
                //boundary_wavs_[p][2*i+1]=boundary_wavs_infact_[2*tempint+1];
                extinfo_[p][2*i] = ((neighbours_[p][2*i] != -1) && !bc[p][2*i]);
                extinfo_[p][2*i+1] = ((neighbours_[p][2*i+1] != -1) && !bc[p][2*i+1]);
            }
        }
        intersection_geometry_centerpatchnumber_.resize(corners_.size(),intpower(9,DIM));
        intersection_geometry_type_.resize(corners_.size(),intpower(9,DIM));
        for (unsigned int i=0; i<DIM; ++i)
        {
            intersection_geometry_orientation_[i].resize(corners_.size(),intpower(9,DIM));
        }
        for (unsigned int i=0; i< corners_.size(); i++)
        {
            for (unsigned int j = 0; j<intpower(9,DIM); j++)
            {
                intersection_geometry_centerpatchnumber_.set_entry(i,j,-1);
            }
        }
    }

    template <class IBASIS, unsigned int DIM>
    QTBasis<IBASIS,DIM>::QTBasis(const Array1D<Point<DIM, int> >& corners,
            const Array1D<FixedArray1D<int,2*DIM> >& neighbours,
            const Array1D<FixedArray1D<int,2*DIM> >& bc)
    : corners_(corners), neighbours_(neighbours), numofbw_((IBASIS::primal_polynomial_degree() + IBASIS::primal_vanishing_moments() -2) / 2)
    {
        assert (corners.size() == neighbours_.size() && neighbours_.size() == bc.size());

        j0_.resize(corners_.size());
        bases_.resize(corners_.size());
        extinfo_.resize(corners_.size());
        //boundary_gens_.resize(corners_.size());
        //boundary_wavs_.resize(corners_.size());
        int tempint;
        typename IBASIS::Index tempind, tempind2;
        IBASIS* b;
        for (unsigned int i=0; i< 4;++i)
        {
            bases_infact_[i]=0;
        }
        for (unsigned int p=0; p < corners_.size();++p)
        {
            for (unsigned int i = 0; i < DIM; i++)
            {
                assert ((bc[p][2*i] < 2) && (bc[p][2*i+1] < 2));
                //cout  << bc[p][2*i] << " " <<  bc[p][2*i+1];
                //tempint = 2 * bc[p][2*i];
                //tempint = bc[p][2*i+1];
                tempint = 2 * bc[p][2*i] + bc[p][2*i+1];
                
                // check whether the corresponding 1d basis already exists
                IBASIS* b = bases_infact_[tempint];
                if (b==0)
                {
                    b = new IBASIS(bc[p][2*i], bc[p][2*i+1]);
                    bases_infact_[tempint] = b;
#if 0
                    // compute number of nonvanishing boundary generators/wavelets
                    if (bc[p][2*i] > 0)
                    {
                        boundary_wavs_infact_[2*tempint] = boundary_gens_infact_[2*tempint] = boundary_wavs_[p][2*i] = boundary_gens_[p][2*i] = 0;
                    }
                    else
                    {
                        tempind = b->first_generator(b->j0());
                        while(true)
                        {
                            if (abs(b->evaluate(0,tempind,0.0)) > 1e-16)
                            {
                                ++boundary_gens_infact_[2*tempint];
                                ++tempind;
                            }
                            else
                            {
                                break;
                            }
                        }
                        assert (boundary_gens_infact_[2*tempint] > 0);
                        tempind = b->first_wavelet(b->j0());
                        while(true)
                        {
                            if (abs(b->evaluate(0,tempind,0.0)) > 1e-16)
                            {
                                ++boundary_wavs_infact_[2*tempint];
                                ++tempind;
                            }
                            else
                            {
                                break;
                            }
                        }
                        assert (boundary_wavs_infact_[2*tempint] > 0);
                    }
                    if (bc[p][2*i+1] > 0)
                    {
                        boundary_wavs_infact_[2*tempint+1] = boundary_gens_infact_[2*tempint+1] = boundary_wavs_[p][2*i+1] = boundary_gens_[p][2*i+1] = 0;
                    }
                    else
                    {
                        tempind = b->first_generator(b->j0());
                        tempind2 = b->last_generator(b->j0());
                        while((abs(b->evaluate(0,tempind,1.0)) < 1e-16) && tempind != tempind2)
                        {
                            ++tempind;
                        }
                        while(true)
                        {
                            ++boundary_gens_infact_[2*tempint+1];
                            if (tempind == tempind2)
                            {
                                break;
                            }
                            else
                            {
                                ++tempind;
                            }
                        }
                        assert (boundary_gens_infact_[2*tempint+1] > 0);

                        tempind = b->first_wavelet(b->j0());
                        tempind2 = b->last_wavelet(b->j0());
                        while((abs(b->evaluate(0,tempind,1.0)) < 1e-16) && tempind != tempind2)
                        {
                            ++tempind;
                        }
                        while(true)
                        {
                            ++boundary_wavs_infact_[2*tempint+1];
                            if (tempind == tempind2)
                            {
                                break;
                            }
                            else
                            {
                                ++tempind;
                            }
                        }
                        assert (boundary_wavs_infact_[2*tempint+1] > 0);
                    }
#endif
                }
                bases_[p][i] = b;
                j0_[p][i] = b->j0();
                //boundary_gens_[p][2*i]=boundary_gens_infact_[2*tempint];
                //boundary_wavs_[p][2*i]=boundary_wavs_infact_[2*tempint];
                //boundary_gens_[p][2*i+1]=boundary_gens_infact_[2*tempint+1];
                //boundary_wavs_[p][2*i+1]=boundary_wavs_infact_[2*tempint+1];
                extinfo_[p][2*i] = ((neighbours_[p][2*i] != -1) && (bc[p][2*i] == 0));
                extinfo_[p][2*i+1] = ((neighbours_[p][2*i+1] != -1) && (bc[p][2*i+1] == 0));
            }
        }
        intersection_geometry_centerpatchnumber_.resize(corners_.size(),intpower(9,DIM));
        intersection_geometry_type_.resize(corners_.size(),intpower(9,DIM));
        for (unsigned int i=0; i<DIM; ++i)
        {
            intersection_geometry_orientation_[i].resize(corners_.size(),intpower(9,DIM));
        }
        for (unsigned int i=0; i< corners_.size(); i++)
        {
            for (unsigned int j = 0; j<intpower(9,DIM); j++)
            {
                intersection_geometry_centerpatchnumber_.set_entry(i,j,-1);
            }
        }
    }

    template <class IBASIS, unsigned int DIM>
    QTBasis<IBASIS,DIM>::~QTBasis()
    {
        for (unsigned int i=0; i<4 ; ++i)
        {
            if ( bases_infact_[i] != 0)
            {
                delete bases_infact_[i];
            }
        }
    }

    template <class IBASIS, unsigned int DIM>
    void 
    QTBasis<IBASIS,DIM>::precompute_firstlast_wavelets()
    {
        // resize level_to_num and num_to_level
        level_to_num_.resize(bases_.size());
        unsigned int temp_int;
        unsigned int num(0); // used as a temporary variable for the moment. Will store the number of the current wavelet later.
        MultiIndex<int,DIM> leveloffset;
        for (unsigned int p = 0; p<bases_.size(); ++p)
        {
            for (unsigned int i = 1 ; i<DIM; ++i)
            {
                leveloffset[i] = 0;
            }
            leveloffset[0] = jmax_-multi_degree(j0_[p]);
            temp_int = leveloffset.number() + 1; // numbering begins at 0
            num += temp_int;
            level_to_num_[p].resize(temp_int);
        }
        num_to_level_.resize(num);
        first_wavelets_.resize(num);
        last_wavelets_.resize(num);
        // determine first level (j_,p_)
        list<int> temp_list;
        for (unsigned int i=0; i< bases_.size(); i++)
        {
            temp_list.push_back(i);
        }
        temp_list.sort(index_cmp< Array1D<MultiIndex<int,DIM> > > (j0_,0));
        int p_= temp_list.front();
        
        // for all patches: store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
        Array1D<FixedArray1D<map<int, int>,DIM> > sizes; 
        sizes.resize(bases_.size());
        for (unsigned int p = 0; p<bases_.size(); ++p)
        {
            for (unsigned int i=0; i < DIM; i++) 
            {
                sizes[p][i][0] = bases_[p][i]->Deltasize(j0_[p][i]); // Number of generators on level j0_[p][i]
                // store all sizes that will be needed up to jmax
                for (unsigned int k = 0; k <= jmax_-multi_degree(j0_[p]); ++k)
                {
                    sizes[p][i][k+1] = bases_[p][i]->Nablasize(j0_[p][i]+k); // Number of Wavelets on level j0+k
                }
            }
        }
        
        // begin computation of first/last wavelets
        // Roadmap:
        // while (levelnorm <= jmax_)
        // {
        //     compute first/last wavelet for current level
        //     increase (j,p), increase levelnorm if needed
        // }
        // code is similar to the operator ++ in qtbasis_index.cpp.
        // Check there for more comments
        MultiIndex<int,DIM> j_(j0_[p_]), e_, k_;
        unsigned int levelnorm(multi_degree(j0_[p_]) ); // norm of the current (minimal) level
        num = 0; // number of current wavelet. Actually the "first wavelet" is a generator 
        unsigned int levelnum(0); //num of the pair (p,j) determined by the ordering of the wavelet indices
        bool done;
        while (levelnorm <= jmax_)
        {
            // compute first/last wavelet for current level
            for (unsigned int i = 0; i < DIM; i++) 
            {
                e_[i] = (j0_[p_][i] != j_[i]);
                k_[i] = (e_[i] == 0) ? (bases_[p_][i]->DeltaLmin()) : (bases_[p_][i]->Nablamin());
            }
            first_wavelets_[levelnum] = Index(j_,e_,k_,p_,num,this);
            temp_int = 1; // number of wavelets on current level
            for (unsigned int i = 0; i<DIM; ++i)
            {
                leveloffset[i] = j_[i]-j0_[p_][i];
                temp_int *= ((e_[i] == 0)? (sizes[p_][i][0] + sizes[p_][i][1] ): (sizes[p_][i][leveloffset[i]+1]) );
            }
            num += temp_int; // this is the number of the next first_wavelet! its 1 too big
            for (unsigned int i = 0; i < DIM; i++) 
            {
                e_[i] = 1;
                k_[i] = (bases_[p_][i]->Nablamax(j_[i]));
            }
            // cout << "levelnum = " << levelnum <<endl;
            // cout << "; first_wavelets_[" << levelnum << "] = " << first_wavelets_[levelnum] <<endl;
            last_wavelets_[levelnum] = Index(j_,e_,k_,p_,num-1,this);
            // cout << "; last_wavelets_[" << levelnum << "] = " << last_wavelets_[levelnum] <<endl;
            num_to_level_[levelnum] = std::make_pair(j_,p_);
            // cout << "; num_to_level_[" << levelnum << "] = ( " << num_to_level_[levelnum].first << ", " << num_to_level_[levelnum].second << endl;
            level_to_num_[p_][leveloffset.number()] = levelnum;
            // cout << "; level_to_num_[" << p_ << "][" << leveloffset.number() << "] = " << level_to_num_[p_][leveloffset.number()] << endl;
            ++levelnum;
            
            // increase (j_,p_)
            // try to increase patch p_
            done = false;
            while ((!done) && (p_ < bases_.size()-1))
            {
                p_++;
                //done = true;
                // check if level j_ exists for the current patch.
                // that is not the same as using > or lex from MultiIndex
                for (unsigned int i = 0; i<DIM; ++i)
                {
                    done = (j_[i] >= j0_[p_][i]);
                    if (!done)
                        break;
                }
            }
            // if needed increase j_
            if (!done)
            {
                if (DIM == 1)
                {
                    ++levelnorm;
                    j_[0]= j_[0]+1;
                    // Assumption: the minimal levels on different patches differ at most by 1:
                    p_ = 0;
                    // Use this code, to allow arbitrary minimal levels:
                    /*
                    for (unsigned int p=0; p < basis_->get_nop(); ++p)
                    {
                        if (j_[0] >= j0_[p][0])
                        {
                            p_ = p;
                        }
                    }
                    */
                }
                else
                {
                    // try to increase the sublevel,
                    // that is increasing j_ without increasing its 1-norm.
                    FixedArray1D<list<int>*, DIM-1> incrementable; // store the patches which allow for an increment
                    Array1D<int> min_suffix_norms; // for each patch store the norm of the suffix (from current to last index)
                    min_suffix_norms.resize(bases_.size());
                    //int min_suffix_norm;

                    //FixedArray1D<list<int>*, DIM-1> prefix;
                    list<int>* temp_list;

                    incrementable[0]=new list<int>;
                    temp_list = new list<int>; // store the patches which allow the prefix up to the current dimension

                    for (unsigned int p=0; p < bases_.size(); p++)
                    {
                        if ( (j0_[p][0] <= j_[0]+1) && (j0_[p][1] < j_[1]) ) // check if decrement at pos=1 and increment at pos =0 are possible
                        {
                            incrementable[0]->push_back(p);
                        }
                    }
                    if (DIM > 2)
                    {
                        for (unsigned int p=0; p < bases_.size(); p++)
                        {
                            if (j0_[p][0] <= j_[0])
                            {
                                temp_list->push_back(p);
                            }
                        }
                    }
                    for (unsigned int i=1; i< DIM-1; i++)
                    {
                        incrementable[i]=new list<int>;
                        for (list<int>::iterator it(temp_list->begin()), itend(temp_list->end()); it != itend;)
                        {
                            if ( (j0_[*it][i] <= j_[i]+1) && (j0_[*it][i+1] < j_[i+1]))
                            {
                                incrementable[i]->push_back(*it);
                            }
                            if (j0_[*it][i] > j_[i])
                            {
                                temp_list->erase(it);
                            }
                            else
                            {
                                ++it;
                            }
                        }
                    }
    
                    // from right to left: search for an index i to decrease j_[i]
                    // done = false;
                    // i == DIM-2
                    if (!incrementable[DIM-2]->empty())
                    {
                        incrementable[DIM-2]->sort();
                        //accept!
                        // only works if minimal differ only by 1
                        p_ = incrementable[DIM-2]->front();
                        j_[DIM-2] += 1;
                        j_[DIM-1] -= 1;

                        done = true;
                    }
                    temp_int = 0;
                    if (!done)
                    {
                        temp_int = j_[DIM-1]; // store j_[i+1] + ... + j_[DIM-1]
                        for (unsigned int p = 0; p < bases_.size();++p)
                        {
                            min_suffix_norms[p]=j0_[p][DIM-1];
                        }
                        for (int i=DIM-3; i>=0; --i)
                        {
                            temp_int += j_[i+1];
                            for (unsigned int p = 0; p < bases_.size();++p)
                            {
                                min_suffix_norms[p]+=j0_[p][i+1];
                            }
                            if (!incrementable[i]->empty())
                            {
                                // we want to increment at position i.
                                // sort increment set according to j0_[p][i+1,...,DIM-2]
                                incrementable[i]->sort(index_cmp_ignoreLastEntry<Array1D<MultiIndex<int,DIM> > > (j0_,i+1));
                                //check min_suffix_norm condition == is the level valid for patch p_ and the current levelnorm
                                for (list<int>::const_iterator it(incrementable[i]->begin()), itend(incrementable[i]->end()); it != itend; ++it)
                                {
                                    // accept if j_[i+1] + ... + j_[DIM-1] > j0_[p][i+1] + ... j0_[p][DIM-1]
                                    if (min_suffix_norms[*it] < temp_int)
                                    {
                                        //accept!
                                        p_ = *it;
                                        j_[i] = j_[i]+1;
                                        for ( unsigned int j = i+1; j <= DIM -2; ++j)
                                        {
                                            j_[j] = j0_[p_][j];
                                        }
                                        j_[DIM-1] = j0_[p_][DIM-1] + temp_int - min_suffix_norms[p_] -1;
                                        done = true;
                                        break;
                                    }
                                }
                            }
                            if (done)
                                break;
                        }
                    }
                    // end of increase-sublevel
                    // "big loop" : increase j_ such that the 1-norm increases by 1
                    if (!done)
                    {
                        assert ((temp_int + j_[0]) == multi_degree(j_));
                        for (unsigned int p = 0; p < bases_.size();++p)
                        {
                            min_suffix_norms[p]+=j0_[p][0];
                            assert (min_suffix_norms[p] == multi_degree(j0_[p]));
                        }
                        temp_int += j_[0]+1; // norm of the new level
                        //temp_int = multi_degree(j_)+1; // norm of the new level
                        temp_list->clear();
                        for (unsigned int p=0; p< bases_.size(); p++)
                        {
                            temp_list->push_back(p);
                        }
                        // We deduce the new patchnumber by sorting the minimal levels
                        temp_list-> sort(index_cmp_ignoreLastEntry< Array1D<MultiIndex<int,DIM> > > (j0_,0));
                        for (list<int>::const_iterator it(temp_list->begin()),itend(temp_list->end());it!=itend;++it)
                        {
                            if (min_suffix_norms[*it] <= temp_int)
                            {
                                //accept this patch
                                j_ = j0_[*it];
                                p_= (*it);
                                j_[DIM-1] = j0_[p_][DIM-1] + temp_int - min_suffix_norms[p_];
                                ++levelnorm;
                                done = true;
                                break;
                            }
                            else
                            {
                                // since we assumed the minimal levels differ by at most 1
                                abort();
                            }
                        }
                    }
                }
            }
            // end of increase (j_,p_)
        }
    }

    
    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::setup_full_collection()
    {
        {
            unsigned int temp_int(0);
            bool jmax_ok = true;
            for (unsigned int p(0);p< get_nop();p++)
            {
                temp_int = std::max(temp_int, multi_degree(j0_[p]) );
                jmax_ok = (jmax_ok && (jmax_ >= temp_int));
            }
            if (!jmax_ok)
            {
                cout << "QTBasis:: setup_full_collection(): specified jmax_ of " << jmax_ << " too low!";
                cout << " Using "<<temp_int << " instead." << endl;
                jmax_ = temp_int;
            }
        }
        full_collection_.resize(last_wavelets_[last_wavelets_.size()-1].number()+1);
        cout << "QTBasis:: total degrees of freedom between j0_ and jmax_ is " << full_collection_.size() << endl;
        cout << "QTBasis:: setting up collection of wavelet indices..." << endl;
        Index ind = first_generator();
        for (int k = 0; k < full_collection_.size(); k++) {
            //cout << "k = " << k << "; ind = " << ind << "; ind.number() = " << ind.number() << endl;
            full_collection_[k] = ind;
            ++ind;
        }
        cout << "done setting up collection of wavelet indices." << endl;
    }



    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::support_local(const Index& lambda, Support& supp) const
    {
        // PERFORMANCE?
        MultiIndex<int,DIM> temp_j=lambda.j(),temp_e=lambda.e(), temp_k=lambda.k();
    	for (unsigned int i(0); i < DIM; i++)
        {
            supp.j[i] = temp_j[i]+temp_e[i];
            //supp.j[i] = lambda.j()[i] + lambda.e()[i];
            bases_[lambda.p()][i]->support(temp_j[i],
                    temp_e[i],
                    temp_k[i],
                    supp.a[i],
                    supp.b[i]);
    	}
    }
    
    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::support_full(const Index& lambda, Support& supp) const
    {
        // PERFORMANCE?
        MultiIndex<int,DIM> temp_j=lambda.j(),temp_e=lambda.e(), temp_k=lambda.k();
    	for (unsigned int i(0); i < DIM; i++)
        {
            supp.j[i] = temp_j[i]+temp_e[i];
            //supp.j[i] = lambda.j()[i] + lambda.e()[i];
            bases_[lambda.p()][i]->support(temp_j[i],
                    temp_e[i],
                    temp_k[i],
                    supp.a[i],
                    supp.b[i]);
            
            if ( (extinfo_[lambda.p()][2*i] == true) 
                    &&
                    ( ( (temp_e[i] == 0) && (temp_k[i] == bases_[lambda.p()][i]->DeltaLmin()) ) 
                  
                      ||
                      ( (temp_e[i] == 1) && (temp_k[i] < numofbw_) ) ) )
            {
                supp.a[i] = - supp.b[i]; // reflect at 0
            }
                    
            if ( (extinfo_[lambda.p()][2*i+1] == true) 
                    &&
                    ( ( (temp_e[i] == 0) && (temp_k[i] == bases_[lambda.p()][i]->DeltaRmax(j0_[lambda.p()][i])) ) 
                      ||
                      ( (temp_e[i] == 1) && (temp_k[i] > (1<<temp_j[i])- 1 - numofbw_) ) ) )
            {
                supp.b[i] = 2 * (1<<temp_j[i]) - supp.a[i]; // reflect at 1
            }
            supp.a[i] += corners_[lambda.p()][i] * (1<< temp_j[i]);
            supp.b[i] += corners_[lambda.p()][i] * (1<< temp_j[i]);
    	}
    }


    template <class IBASIS, unsigned int DIM>
    double
    QTBasis<IBASIS,DIM>::evaluate_check(const unsigned int derivative,
                                      const int lambdanum,
                                      const Point<DIM> x) const
    {
        double value = 1.0;
        Index lambda (full_collection_[lambdanum]);
        int p = lambda.p(), temp_e, temp_j, temp_k;
        double relative;
        for (unsigned int i = 0; i < DIM; i++) // loop through components of the tensor product
        {
            relative = x[i]-corners_[p][i];
            if (relative < 0)
            {
                temp_e = lambda.e()[i];
                temp_k = lambda.k()[i];
                // check whether the basis function is extended in the current direction:
                
                if ((extinfo_[lambda.p()][2*i] == true)  
                    &&
                   (((temp_e == 0) && (temp_k == bases_[p][i]->DeltaLmin()) ) || ((temp_e == 1) && (temp_k < numofbw_ ))))
                {
                    value *= bases_[p][i]->evaluate(derivative,
                                        lambda.j()[i], temp_e, temp_k,
                                         - relative); // reflection at 0
                }
                else
                {
                    value = 0;
                    break;
                }
            }
            else if (relative > 1)
            {
                temp_e = lambda.e()[i];
                temp_j = lambda.j()[i];
                temp_k = lambda.k()[i];
                // check whether the basis function is extended in the current direction:
                
                if ((extinfo_[lambda.p()][2*i+1] == true) 
                    &&
                    (((temp_e == 0) && (temp_k == bases_[p][i]->DeltaRmax(j0_[p][i]))) || ((temp_e == 1) && (temp_k > (1<<temp_j)- 1-numofbw_))))
                {
                    value *= bases_[p][i]->evaluate(derivative,
                                        temp_j, temp_e, temp_k,
                                        2 - relative);
                }
                else
                {
                    value = 0;
                    break;
                }
            }
            else
            {
                value *= bases_[p][i]->evaluate(derivative,
                                             lambda.j()[i], lambda.e()[i], lambda.k()[i],
                                             relative);
            }
        }
        return value;
    }

    template <class IBASIS, unsigned int DIM>
    double
    QTBasis<IBASIS,DIM>::evaluate_simple(const unsigned int derivative,
                                             const int lambdanum,
                                             const Point<DIM> x) const
    {
        double value = 1.0, relative;
        Index lambda (full_collection_[lambdanum]);
        int p = lambda.p();
        for (unsigned int i = 0; i < DIM; i++) // loop through components of the tensor product
        {
            relative = x[i]-corners_[p][i];
            value *= bases_[p][i]->evaluate(derivative,
                                            lambda.j()[i], lambda.e()[i], lambda.k()[i],
                                            (relative < 0) ? (-relative) : ((relative > 1)? (2-relative): relative) );
        }
        return value;
    }
    
    template <class IBASIS, unsigned int DIM>
    double
    QTBasis<IBASIS,DIM>::integrate(const Function<DIM>* f,
                                   const unsigned int lambdanum) const
    {        
        // f(v) = \int_0^1 g(t)v(t) dt
        double r = 0;
        // first compute supp(psi_lambda)
        Support supp_loc;
// CLEANUP !!        
       //,supp_abs;
        Index lambda(full_collection_[lambdanum]);
        support_local(lambda, supp_loc);
        FixedArray1D<bool,DIM> extleft,extright; // psi_lambda is extended left or right in the current direction?
        MultiIndex<int,DIM> temp_j=lambda.j(),temp_e=lambda.e(), temp_k=lambda.k();
        int temp_p(lambda.p());

    	for (unsigned int i(0); i < DIM; i++)
        {
// CLEANUP
            extleft[i] = false;
            extright[i] = false;
            //supp_abs.j = supp_loc.j;
            //supp_abs.a = supp_loc.a;
            //supp_abs.b = supp_loc.b;
            if ( (extinfo_[temp_p][2*i] == true) 
                    &&
                    ( ( (temp_e[i] == 0) && (temp_k[i] == bases_[temp_p][i]->DeltaLmin()) ) 
                  
                      ||
                      ( (temp_e[i] == 1) && (temp_k[i] < numofbw_) ) ) )
            {
                extleft[i] = true;
// CLEANUP
                assert (supp_loc.a[i] == 0);
                //supp_abs.a[i] = - supp_abs.b[i]; // reflect at 0
            }
                    
            if ( (extinfo_[temp_p][2*i+1] == true) 
                    &&
                    ( ( (temp_e[i] == 0) && (temp_k[i] == bases_[temp_p][i]->DeltaRmax(j0_[temp_p][i])) ) 
                      ||
                      ( (temp_e[i] == 1) && (temp_k[i] > (1<<temp_j[i])- 1 - numofbw_) ) ) )
            {
                extright[i] = true;
//CLEANUP
                assert (supp_loc.b[i] == (1<<(temp_j[i]+temp_e[i]) ));
                //supp_abs.b[i] = 2 * (1<<temp_j[i]) - supp_abs.a[i]; // reflect at 1
            }
            //supp_abs.a[i] += corners_[lambda.p()][i] * (1<< temp_j[i]);
            //supp_abs.b[i] += corners_[lambda.p()][i] * (1<< temp_j[i]);
        }
        
        // setup Gauss points and weights for a composite quadrature formula:
        const int N_Gauss = 5;
        FixedArray1D<double,DIM> h;
        for (unsigned int i=0; i<DIM;i++)
        {
            h[i]=ldexp(1.0, -supp_loc.j[i]); // granularity for the quadrature
        }
        //FixedArray1D<Array1D<double>,DIM> gauss_points_abs; // absolute = real support,
        FixedArray1D<Array1D<double>,DIM> gauss_points, gauss_weights, v_values; //  local = without shift or extension, i.e., on the unit cube
        for (unsigned int i = 0; i < DIM; i++) 
        {
            gauss_points[i].resize(N_Gauss*(supp_loc.b[i]-supp_loc.a[i]));
            gauss_weights[i].resize(N_Gauss);
            /*
            gauss_points_abs[i].resize(N_Gauss*(supp_abs.b[i]-supp_abs.a[i]));
            for (int patch = supp_abs.a[i]; patch < supp_abs.b[i]; patch++)
            {
                for (int n = 0; n < N_Gauss; n++) 
                {
                    gauss_points_abs[i][(patch-supp_abs.a[i])*N_Gauss+n]
                            = h[i]*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                }
            }
            */
            for (int patch = supp_loc.a[i]; patch < supp_loc.b[i]; patch++)
            {
                for (int n = 0; n < N_Gauss; n++) 
                {
                    gauss_points[i][(patch-supp_loc.a[i])*N_Gauss+n]
                            = h[i]*(2*patch+1+GaussPoints[N_Gauss-1][n])/2.;
                }
            }
            for (int n = 0; n < N_Gauss; n++) 
            {
                gauss_weights[i][n]
                        = h[i]*GaussWeights[N_Gauss-1][n];
            }
        }
        // compute the point values of the integrand
        // - we use that it is a tensor product
        // - we only compute values on the unit cube
        for (unsigned int i = 0; i < DIM; i++)
            bases_[temp_p][i]->evaluate(0,
                                 temp_j[i],temp_e[i],temp_k[i],
                                 gauss_points[i],
                                 v_values[i]);
        // iterate over all points and sum up the integral shares
        // we iterate over each point only once, no matter how psi_\lambda is extended
        // thus we have to introduce an iterator over the relevant patches
        // attention: do not count the boundary point twice
        // iterator for the points on the unit cube:
        int index[DIM]; 
        // iterator for the patches lambda is extended to:
        // since we assume that the one dimensional wavelets are only extended to at most one side, an bool array is sufficient
        bool cube_index[DIM]; 
        for (unsigned int i = 0; i < DIM; i++)
        {
                index[i] = 0; // the first point in this direction
                cube_index[i] = false; // the first extension direction ("left" or "middle")
        }
        Point<DIM> x_loc;
        double share(0);
        bool exit = false;
        while (!exit) 
        {
            share = 0;
            while (true) // collect function values f(x) on all cubes lambda is extended to
            {
                for (unsigned int i = 0; i < DIM; i++)
                {
                    // extleft[i] == false, cube_index[i] == false => Mitte
                    // extleft[i] == false, cube_index[i] == true  => Rechts
                    // extleft[i] == true,  cube_index[i] == false => Links
                    // extleft[i] == true,  cube_index[i] == true  => Mitte
                    if (cube_index[i] == extleft[i])
                    { // Mitte
                        x_loc[i] = corners_[temp_p][i] + gauss_points[i][index[i]];
                    }
                    else
                    {
                        if ((cube_index[i] == false) && (extleft[i] == true))
                        {
                            // Links
                            x_loc[i] = corners_[temp_p][i] - gauss_points[i][index[i]];
                        }
                        else
                        {
                            // Rechts
                            x_loc[i] = corners_[temp_p][i] + 2. - gauss_points[i][index[i]];
// CLEANUP                  
                            /*
                            cout << " in qtbasis.cpp :: integrate: extension to the right! " <<
                                    "gauss_points[" << i<< "][" << index[i] << "] = " << gauss_points[i][index[i]] <<
                                    "; x_loc[" << i << "] = " << x_loc[i] << endl;
                             */
                        }
                    }
                }
                share += f->value(x_loc);
                // "++cube_index"
                //exit = false;
                for (unsigned int i = 0; i < DIM; i++) 
                {
// CLEANUP                    
                    /*
                    if ( ( (extleft[i] == false) && (extright[i] == false) )
                        ||
                            (( (extleft[i] == true) || (extright[i] == true) ) && (cube_index[i] == true)) )
                    {
                        cube_index[i] = false;
                        exit = (i == DIM-1);
                    }
                    else */
                    
                    if ((cube_index[i] == false) && ( (extleft[i] == true) || (extright[i] == true) ) )
                    {
                        cube_index[i] = true;
                        break;
                    }
                    else
                    {
                        cube_index[i] = false;
                        exit = (i == DIM-1);
                    }
                }
                if (exit) break;
            }
            
            for (unsigned int i = 0; i < DIM; i++)
            {
                share *= gauss_weights[i][index[i]%N_Gauss] * v_values[i][index[i]];
                
            }
            r += share;
            // "++index"
            //exit = false;
            for (unsigned int i = 0; i < DIM; i++) 
            {
                if (index[i] == N_Gauss*(supp_loc.b[i]-supp_loc.a[i])-1) 
                {
                    index[i] = 0;
                    exit = (i == DIM-1);
                } else 
                {
                    index[i]++;
                    exit = false;
                    break;
                }
            }
            //if (exit) break
        }
        return r;
    }
    
    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::expand(const Function<DIM>* f,
                                const bool primal,
                                //const MultiIndex<int,DIM> jmax,
                                InfiniteVector<double,Index>& coeffs) const
    {
        assert(primal == false); // only integrate against primal wavelets and generators
        for (unsigned int k=0; k< degrees_of_freedom(); ++k)
        //for (Index lambda = first_generator(), lambda_end(last_wavelet(jmax));;++lambda)
        {
            const double coeff = integrate(f, full_collection_[k]);
            //const double coeff = integrate(f, lambda);
            if (fabs(coeff)>1e-15)
            {
                coeffs.set_coefficient(full_collection_[k], coeff);
                //coeffs.set_coefficient(lambda, coeff);
            }
        }
    }

    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::expand(const Function<DIM>* f,
                                const bool primal,
                                InfiniteVector<double,int>& coeffs) const
    {
        assert(primal == false); // only integrate against primal wavelets and generators
        for (unsigned int k=0; k< degrees_of_freedom(); ++k)
        //for (Index lambda = first_generator(), lambda_end(last_wavelet(jmax));;++lambda)
        {
            const double coeff = integrate(f, k);
            //const double coeff = integrate(f, lambda);
            if (fabs(coeff)>1e-15)
            {
                coeffs.set_coefficient(k, coeff);
                //coeffs.set_coefficient(lambda, coeff);
            }
        }
    }

    template <class IBASIS, unsigned int DIM>
    typename QTBasis<IBASIS,DIM>::Index
    QTBasis<IBASIS,DIM>::first_generator() const
    {
        list<int> temp_list;
        for (unsigned int i=0; i< get_nop(); i++)
        {
            temp_list.push_back(i);
        }
        temp_list.sort(index_cmp< Array1D<MultiIndex<int,DIM> > > (j0_,0));

        int p= temp_list.front();
        typename Index::type_type e;
        typename Index::translation_type k;
        for (unsigned int i = 0; i < DIM; i++) {
            e[i]=0;
            k[i]=bases_[p][i]->DeltaLmin();
        }
        return Index(j0_[p], e, k, p, 0, this);
    }
    
    
    
    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::get_cached_onedim_intersections(const unsigned int lambda_basisnum,
                const unsigned int lambda_j,
                const unsigned int lambda_e,
                const unsigned int lambda_k,
                const bool reflected,
                const unsigned int mu_basisnum,
                const unsigned int mu_j,
                const bool generators,
                int& kmin,
                int& kmax)
    {
        //typename IBASIS::Index lam1(3,0,0,bases_infact_[0]);
        //typename IBASIS::Index lam2(3,1,0,bases_infact_[0]);
        //cout << lam1.number() << endl;
        //cout << lam2.number() << endl;
        typename IBASIS::Index lambda(lambda_j,lambda_e,lambda_k,bases_infact_[lambda_basisnum]);
// cleanup
        int lam_k1, lam_k2,temp_i1, temp_i2;
        bases_infact_[lambda_basisnum]->support(lambda_j,lambda_e,lambda_k,lam_k1,lam_k2);
        if (reflected)
        {
            temp_i1 = 1<< (lambda_j + lambda_e);
            temp_i2 = temp_i1-lam_k2;
            lam_k2 = temp_i1-lam_k1;
            lam_k1 = temp_i2;
        }
        get_intersecting_wavelets_on_level(*bases_infact_[mu_basisnum],
                        lambda_j,
                        lambda_e,
                        lam_k1,
                        lam_k2,
                        mu_j, 
                        generators, 
                        temp_i1, 
                        temp_i2);
        
        if (reflected)
        {
            if (generators)
            {
                if (reflected_gens[lambda_basisnum*4+mu_basisnum].count(lambda.number()) == 1)
                {
                    kmin = reflected_gens[lambda_basisnum*4+mu_basisnum][lambda.number()].first;
                    kmax = reflected_gens[lambda_basisnum*4+mu_basisnum][lambda.number()].second;
                    assert (kmin == temp_i1);
                    assert (kmax == temp_i2);
                    return;
                }
            }
            else
            {
                if (reflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()].count(mu_j-bases_infact_[mu_basisnum]->j0()) == 1 )
                {
                    kmin = reflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()][mu_j-bases_infact_[mu_basisnum]->j0()].first;
                    kmax = reflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()][mu_j-bases_infact_[mu_basisnum]->j0()].second;
                    assert (kmin == temp_i1);
                    assert (kmax == temp_i2);
                    return;
                }
            }
        }
        else
        {
            if (generators)
            {
                if (nonreflected_gens[lambda_basisnum*4+mu_basisnum].count(lambda.number()) == 1)
                {
                    //cout << lambda.number() << endl;
                    kmin = nonreflected_gens[lambda_basisnum*4+mu_basisnum][lambda.number()].first;
                    kmax = nonreflected_gens[lambda_basisnum*4+mu_basisnum][lambda.number()].second;
                    assert (kmin == temp_i1);
                    assert (kmax == temp_i2);
                    return;
                }
            }
            else
            {
                if (nonreflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()].count(mu_j-bases_infact_[mu_basisnum]->j0()) == 1 )
                {
                    //cout << "lambda_basisnum*4+mu_basisnum = " << lambda_basisnum*4+mu_basisnum << endl;
                    //cout << "lambda = " << lambda << "; lambda.number() = " << lambda.number() << endl;
                    //cout << "mu_j-bases_infact_[mu_basisnum]->j0() = " << mu_j << " - " << bases_infact_[mu_basisnum]->j0() << " = " << mu_j-bases_infact_[mu_basisnum]->j0() << endl;
                    kmin = nonreflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()][mu_j-bases_infact_[mu_basisnum]->j0()].first;
                    kmax = nonreflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()][mu_j-bases_infact_[mu_basisnum]->j0()].second;
                    assert (kmin == temp_i1);
                    assert (kmax == temp_i2);
                    //cout << "kmin = " << kmin << endl;
                    //cout << "kmax = " << kmax << endl;
                    return;
                }
            }
        }
        // value is not cached. compute it!
        int temp_k1,temp_k2;
        bases_infact_[lambda_basisnum]->support(lambda_j,lambda_e,lambda_k,temp_k1,temp_k2);
        if (reflected)
        {
            int temp_i3(1<< (lambda_j + lambda_e)), temp_i4;
            temp_i4 = temp_i3-temp_k2;
            temp_k2 = temp_i3-temp_k1;
            temp_k1 = temp_i4;
        }
        get_intersecting_wavelets_on_level(*(bases_infact_[mu_basisnum]),
                lambda_j,
                lambda_e,
                temp_k1,
                temp_k2,
                mu_j, generators, kmin,kmax);
// cleanup        
        assert (kmin == temp_i1);
        assert (kmax == temp_i2);
        
        if (reflected)
        {
            if (generators)
            {
                reflected_gens[lambda_basisnum*4+mu_basisnum][lambda.number()] = std::pair<int,int> (kmin,kmax);
            }
            else
            {
                reflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()][mu_j-bases_infact_[mu_basisnum]->j0()] = std::pair<int,int> (kmin,kmax);
            }
        }
        else
        {
            if (generators)
            {
                nonreflected_gens[lambda_basisnum*4+mu_basisnum][lambda.number()] = std::pair<int,int> (kmin,kmax);
            }
            else
            {
                nonreflected_wavs[lambda_basisnum*4+mu_basisnum][lambda.number()][mu_j-bases_infact_[mu_basisnum]->j0()] = std::pair<int,int> (kmin,kmax);
            }
        }
//cleanup
        get_cached_onedim_intersections(lambda_basisnum,
                lambda_j,
                lambda_e,
                lambda_k,
                reflected,
                mu_basisnum,
                mu_j,
                generators,
                kmin,
                kmax);
        assert (kmin == temp_i1);
        assert (kmax == temp_i2);
    }
    
    template <class IBASIS, unsigned int DIM>
    bool
    QTBasis<IBASIS,DIM>::get_LMR_info(const unsigned int lambdanum, 
                //const unsigned int levelnum,
                //const MultiIndex<int,DIM> mu_j,
                const unsigned int mu_p,
                MultiIndex<unsigned int, DIM>& intinfo) const
    {
        // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
        Index lambda(get_wavelet(lambdanum));
        unsigned int lambda_p(lambda.p());
        MultiIndex<int,DIM> lambda_j(lambda.j()), lambda_e(lambda.e()), lambda_k(lambda.k());
        //compute intinfo. Similar code as in "integrate"
        int temp_i1, temp_i2;
        for (unsigned int i(0); i < DIM; i++)
        {
            temp_i1 = 1;
            if ( (extinfo_[lambda_p][2*i] == true)  // patch von lambda, linker Rand
                    &&
                    ( ( (lambda_e[i] == 0) && (lambda_k[i] == bases_[lambda_p][i]->DeltaLmin()) )  // Annahme: nur der erste und letzte Generator sind am Rand nichttrivial und müssen fortgesetzt werden
                      ||
                      ( (lambda_e[i] == 1) && (lambda_k[i] < numofbw_) ) ) ) // Annahme: first_wavelet.k  == 0
            {
                temp_i1 = 0; // lambda wird in dieser dimension links fortgesetzt
            }
            if ( (extinfo_[lambda_p][2*i+1] == true) // patch von lambda, rechter Rand
                    &&
                    ( ( (lambda_e[i] == 0) && (lambda_k[i] == bases_[lambda_p][i]->DeltaRmax(j0_[lambda_p][i])) ) 
                      ||
                      ( (lambda_e[i] == 1) && (lambda_k[i] > (1<<lambda_j[i])- 1 - numofbw_) ) ) )
            {
// CLEANUP
                assert (temp_i1 == 1); //otherwise lambda would be extended left and right
                temp_i1 = 2; // lambda muss rechts fortgesetzt werden
            }
            temp_i2 = corners_[lambda_p][i] - corners_[mu_p][i];
            if (abs(temp_i2) > 1)
            {
                // no intersection possible
                return false;
            }
            if (temp_i2 == 0)
            {
                if ( (temp_i1==0) && (extinfo_[mu_p][2*i]) )
                {
                    // LL
                    intinfo[i] = 0;
                }
                else
                {
                    if ( (temp_i1==2) && (extinfo_[mu_p][2*i+1]) )
                    {
                        // RR
                        intinfo[i] = 8;
                    }
                    else
                    {
                        // MM
                        intinfo[i] = 4;
                    }
                }
            }
            else
            {
                if (temp_i2 > 0) // mu_p is left of lambda_p
                {
                    if (temp_i1 == 0)
                    {
                        // (extinfo[mup][2*i+1] == true) ? LR : LM
                        intinfo[i] = (extinfo_[mu_p][2*i+1] == true) ? 2 : 1;
                    }
                    else
                    {
                        if ( (temp_i1 == 1) && (extinfo_[mu_p][2*i+1] == true) )
                        {
                            // MR
                            intinfo[i] = 5;
                        }
                        else
                        {
                            //no intersection
                            return false;
                        }
                    }
                }
                else
                {
                    if (temp_i1 == 2)
                    {
                        // (extinfo[mup][2*i] == true) ? RL : RM
                        intinfo[i] = (extinfo_[mu_p][2*i] == true) ? 6 : 7;
                    }
                    else
                    {
                        if ( (temp_i1 == 1) && (extinfo_[mu_p][2*i] == true) )
                        {
                            // ML
                            intinfo[i] = 3;
                        }
                        else
                        {
                            // no intersection
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }
    
     
    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::get_intersection_geometry(const unsigned int lambda_p, 
                const MultiIndex<unsigned int, DIM> intinfo, 
                unsigned int& type, 
                int& centerpatchnumber, 
                MultiIndex<bool, DIM> & orientation) 
    {
        unsigned int intinfonumber(0);
        for (unsigned int i=0; i<DIM; ++i)
        {
            intinfonumber *= 9;
            intinfonumber += intinfo[i];
        }
        
        if (intersection_geometry_centerpatchnumber_.get_entry(lambda_p, intinfonumber) != -1)
        {
            centerpatchnumber = intersection_geometry_centerpatchnumber_.get_entry(lambda_p, intinfonumber);
            type = intersection_geometry_type_.get_entry(lambda_p, intinfonumber);
            for (unsigned int i=0; i<DIM; ++i)
            {
                orientation[i] = intersection_geometry_orientation_[i].get_entry(lambda_p, intinfonumber);
            }
        }
        else
        {
            /*
             * intinfo[i] orientation centerpatch
             * 0 LL refl left_of_lambda
             * 1 LM ok left_of_lambda
             * 2 LR ok left_of_lambda
             * 3 ML refl same_as_lambda
             * 4 MM ok same_as_lambda
             * 5 MR refl same_as_lambda
             * 6 RL refl same_as_lambda
             * 7 RM ok right_of_lambda
             * 8 RR ok same_as_lambda
             */
            
            type = 0;
            int temp_i1;
            centerpatchnumber = lambda_p;
            for (int i=DIM-1; i >=0; --i)
            {
                temp_i1 = intinfo[i];
                type *= 2;
                type += ((temp_i1%2 == 0) && (temp_i1 !=4)) ? 1:0 ;
                orientation[i] = (temp_i1 != 0) && ((temp_i1 != 3) && ((temp_i1 != 5) && (temp_i1 != 6))); // not LL, ML,MR, RL
                centerpatchnumber = (temp_i1 < 3)? (neighbours_[centerpatchnumber][2*i]) // look left in the current dimension
                                                 : ( (temp_i1 == 7)? (neighbours_[centerpatchnumber][2*i+1]) : (centerpatchnumber) ); // RM? then look right!
            }
            /*
             * 
            int temp_i1(intinfo[0]);
            int temp_i2(intinfo[1]);
            // intinfo[0] == LL, LR, RL, RR == 0,2,6,8
            bool temp_b1 ( ((temp_i1%2 == 0) && (temp_i1 !=4)) );
            bool temp_b2 ( ((temp_i2%2 == 0) && (temp_i2 !=4)) );
            // if(temp_b1 == true) -> geometry 1,3,4,5,6,7 else 0,2
            // if(temp_b2 == true) -> geometry 1,3,4,5,6,7 else 0,2
            type = temp_b1? (temp_b2? (7)
                                    : (1) )
                          :( temp_b2? (2)
                                    : (0) );
            // if (type == 7) ... check 3-6
            // orientation[i] == true if mu and centerpatch are on the same line (of the quadrangulation), i.e.,
            // intinfo[i] == !LL, LM, LR, !ML, MM, 
                    
            orientation[0] = (temp_i1 != 0) && ((temp_i1 != 3) && ((temp_i1 != 5) && (temp_i1 != 6))); // not LL, ML,MR, RL
            orientation[1] = (temp_i2 != 0) && ((temp_i2 != 3) && ((temp_i2 != 5) && (temp_i2 != 6))); // not LL, ML,MR, RL
            centerpatchnumber = (temp_i1 < 3)? (neighbours_[lambda_p][0])
                                             : ( (temp_i1 == 7)? (neighbours_[lambda_p][1]) : (lambda_p) );
            centerpatchnumber = (temp_i2 < 3)? (neighbours_[centerpatchnumber][2]) 
                                             : ( (temp_i2 == 7)? (neighbours_[centerpatchnumber][3]) : (centerpatchnumber) );
            */
            intersection_geometry_centerpatchnumber_.set_entry(lambda_p, intinfonumber , centerpatchnumber);
            intersection_geometry_type_.set_entry(lambda_p, intinfonumber, type);
            for (unsigned int i=0; i<DIM; ++i)
            {
                intersection_geometry_orientation_[i].set_entry(lambda_p, intinfonumber, orientation[i]);
            }
        }
    }
    
    template <class IBASIS, unsigned int DIM>
    void
    QTBasis<IBASIS,DIM>::get_onedim_intersections(const unsigned int intinfo_i,
            const unsigned int lami_j, 
            const unsigned int lami_e, 
            const int lami_k,
            const unsigned int lami_basisnum,
            const bool min_type_i,
            const unsigned int mui_j,
            const unsigned int mui_basisnum,
            int& kmingen_i,
            int& kmaxgen_i,
            int& kminwav_i,
            int& kmaxwav_i,
            bool& gen_intersection,
            bool& wav_intersection)
    {
        // it is possible that there exists a nontrivial set of intersecting mus
        // however, if lambda[i] is "in the middle of the interval" this is not guaranteed.
        // Example:
        // mu[i] is extended, i.e., intinfo[i]=MR, then it is possible that none of the extended boundary wavelets mu[i] can reach lambda[i]
        // analogous if mu[i] is extended to the left, ...
        // support of boundary generators is smaller than of boundary wavelets (on the same level). 
        // Therefore it may happen that there exists an intersecting boundary wavelet, but not an intersecting boundary generator.

        // compute kmingen,kmaxgen,kminwav,kmaxwav
        
        // reflect support of lambda if needed, i.e.,
        // transform supp_lam into support w.r.t. the segment of mu
        // {LL, LM, LR, ML, MM, MR, RL, RM, RR} = {0,1,...,8}
        
                
        //min_type=(mui_j==j0_[mu_p][i]) ? true:false; // remember whether we are on the minimal levels to avoid calls of basis.j0()
        
        ////lambdabasisnum = ((bc[lambda_p][2*i])?0:2) + ((bc[lambda_p][2*i+1])?0:1);
        ////mubasisnum = ((bc[mu_p][2*i])?0:2) + ((bc[mu_p][2*i+1])?0:1);
        gen_intersection = false;
        wav_intersection = false;
        const int scale_lami(lami_j + lami_e);
        const int scale_mui_wav(mui_j + 1);
        int max_scale;
        int temp_i1, temp_i2;
        int lami_k1,lami_k2,mui_k1,mui_k2;
        
        if ((intinfo_i == 3) || (intinfo_i == 5))
        {
            bases_infact_[lami_basisnum]->support(lami_j,
                            lami_e,
                            lami_k,
                            lami_k1,
                            lami_k2);
            // support lami = 2^{-(j+e)}[lami_k1,lami_k2]
            temp_i1 = (1 << scale_lami)-lami_k2;
            lami_k2 = (1 << scale_lami)-lami_k1;
            lami_k1 = temp_i1;
        }
        
        if (min_type_i == true)
        {
            const int scale_mui_gen = mui_j; // generator case
            switch (intinfo_i) 
            {
                case 3:
                case 5:
                    gen_intersection = false;
                    if (intinfo_i == 3) // ML
                    {
                        bases_infact_[mui_basisnum]->support(bases_infact_[mui_basisnum]->first_generator(mui_j),
                            mui_k1,
                            mui_k2);
                    }
                    else // MR
                    {
                        bases_infact_[mui_basisnum]->support(bases_infact_[mui_basisnum]->last_generator(mui_j),
                            mui_k1,
                            mui_k2);
                    }
                    //const int scale_lambda = lambda.j()[i] + lambda.e()[i];

                    max_scale = std::max(scale_lami,scale_mui_gen);

                    temp_i1 = mui_k2 <<(max_scale-scale_mui_gen);
                    temp_i2 = lami_k1 <<(max_scale-scale_lami);
    // clean up                            
                    //if (temp_i2 >= temp_i1)
                    //{
                        // no intersection with generators
                        //assert(!temp_b); // intersection with wavelets but not with a generator? that is possible! Support of boundary gens is smaller than of boundary wavs
                        //temp_b = (temp_b || false);
                    //}
                    if (temp_i2 < temp_i1)    
                    {
                        temp_i1 = mui_k1 <<(max_scale-scale_mui_gen);
                        temp_i2 = lami_k2 <<(max_scale-scale_lami);
    // clean up                            
                        //if (temp_i2 <= temp_i1)
                        //{
                            // no intersection with generators
                            //assert(!temp_b); // intersection with wavelets but not with a generator? that is possible! Support of boundary gens is smaller than of boundary wavs
                            //temp_b = (temp_b || false);
                        //}
                        if (temp_i2 > temp_i1)
                        {
                            // reflected lambda intersects with boundary generator
                            // only one generator is nontrivial at the boundary (and thus extended)
                            if (intinfo_i == 3) // ML
                            {
                                kmingen_i = bases_infact_[mui_basisnum]->DeltaLmin();
                                kmaxgen_i = kmingen_i;
    // clean up
                            assert(bases_infact_[mui_basisnum]->first_generator(mui_j).k() == bases_infact_[mui_basisnum]->DeltaLmin());
                            }
                            else // MR
                            {
                                kmingen_i = bases_infact_[mui_basisnum]->DeltaRmax(mui_j);
                                kmaxgen_i = kmingen_i;
    // clean up
                                assert(bases_infact_[mui_basisnum]->last_generator(mui_j).k() == bases_infact_[mui_basisnum]->DeltaRmax(mui_j));
                            }
    // clean up                            
                            //assert(temp_b); // no intersection with any wavelet but with a generator? unlikely!
                            gen_intersection = true;
                        }
                    }
                    break;
                default:
                    get_cached_onedim_intersections (
                        lami_basisnum,
                        lami_j,
                        lami_e,
                        lami_k,
                        !(((intinfo_i == 0) || (intinfo_i == 4) ) || (intinfo_i == 8) ),
                        mui_basisnum,
                        mui_j,
                        true,
                        kmingen_i,
                        kmaxgen_i);
                    gen_intersection = true;
                    break;
            }
        }
        
        switch (intinfo_i) {
            case 3:
            case 5:
                wav_intersection = false;
                // 3: ML; 5: MR 
                
                if (intinfo_i == 3)
                {
                    bases_infact_[mui_basisnum]->support(bases_infact_[mui_basisnum]->first_wavelet(mui_j),
                            mui_k1,
                            mui_k2);
                }
                else
                {
                    bases_infact_[mui_basisnum]->support(bases_infact_[mui_basisnum]->last_wavelet(mui_j),
                            mui_k1,
                            mui_k2);
                }
                // reflect supp lambda and check whether it intersects with a left boundary wavelet
                // similar to p_basis:: intersecting_supports:

                //int scale_lami = lami_j + lami_e;
                //int scale_mui_wav = mui_j + 1; // wavelet case

                

                max_scale = std::max(scale_lami,scale_mui_wav);

                temp_i1 = mui_k2 <<(max_scale-scale_mui_wav);
                temp_i2 = lami_k1 <<(max_scale-scale_lami);
                
                if (temp_i2 < temp_i1)
                {
// CLEANUP
                    temp_i1 = mui_k1 <<(max_scale-scale_mui_wav);
                    temp_i2 = lami_k2 <<(max_scale-scale_lami);
                    if (temp_i2 > temp_i1)
                    {
                        // reflected lambda intersects with boundary wavelet
                        // all boundary wavelets have the same support!
                        if (intinfo_i == 3) // ML
                        {
                            kminwav_i = 0; //
                            kmaxwav_i = numofbw_-1;
// clean up
                        assert(bases_infact_[mui_basisnum]->first_wavelet(mui_j).k() == 0);
                        }
                        else // MR
                        {
                            kmaxwav_i = (1<<mui_j)-1;
                            kminwav_i = kmaxwav_i -numofbw_+1;
// clean up
                            assert(bases_infact_[mui_basisnum]->last_wavelet(mui_j).k() == kmaxwav_i);
                        }
                        wav_intersection = true;
                    }
                }
                break;
            default:
                wav_intersection = true;
                get_cached_onedim_intersections ( lami_basisnum,
                    lami_j,
                    lami_e,
                    lami_k,
                    !(((intinfo_i == 0) || (intinfo_i == 4) ) || (intinfo_i == 8) ),
                    mui_basisnum,
                    mui_j,
                    false,
                    kminwav_i,
                    kmaxwav_i);
// CLEANUP                
                int kminwav_i2, kmaxwav_i2;
                get_cached_onedim_intersections ( lami_basisnum,
                    lami_j,
                    lami_e,
                    lami_k,
                    !(((intinfo_i == 0) || (intinfo_i == 4) ) || (intinfo_i == 8) ),
                    mui_basisnum,
                    mui_j,
                    false,
                    kminwav_i2,
                    kmaxwav_i2);
                assert (kminwav_i == kminwav_i2);
                assert (kmaxwav_i == kmaxwav_i2);
                
                break;
        }
// CLEANUP
        if (min_type_i)
        {
            assert (!gen_intersection || wav_intersection);
        }
    }
            
    template <class IBASIS, unsigned int DIM>
    bool
    QTBasis<IBASIS,DIM>::intersecting_wavelets(const unsigned int lambdanum, 
                const unsigned int levelnum,
                MultiIndex<int, DIM>& kmingen,
                MultiIndex<int, DIM>& kmaxgen,
                MultiIndex<int, DIM>& kminwav,
                MultiIndex<int, DIM>& kmaxwav,
                MultiIndex<bool, DIM>& gen_intersection,
                MultiIndex<bool, DIM>& wav_intersection,
                MultiIndex<unsigned int, DIM>& intinfo)
    {
        unsigned int mu_p;
        MultiIndex<int,DIM> mu_j;
        get_level(levelnum,mu_j,mu_p);
        bool temp_b = get_LMR_info(lambdanum, 
                mu_j,
                mu_p,
                intinfo);
        if (!temp_b) 
            return temp_b;
        
        Index lambda(get_wavelet(lambdanum));
        MultiIndex<int,DIM> lambda_j(lambda.j()), lambda_e(lambda.e()), lambda_k(lambda.k());
        int lambda_p(lambda.p());
        bool min_type_i;
        for (unsigned int i=0; i< DIM; i++)
        {
            min_type_i = (mu_j[i]==j0_[mu_p][i]); // remember whether we are on the minimal levels to avoid calls of basis.j0()
            get_onedim_intersections(intinfo[i],
                    lambda_j[i],
                    lambda_e[i],
                    lambda_k[i],
                    (((bc_[lambda_p][2*i])?0:2) + ((bc_[lambda_p][2*i+1])?0:1)),
                    min_type_i,
                    mu_j[i],
                    (((bc_[mu_p][2*i])?0:2) + ((bc_[mu_p][2*i+1])?0:1)),
                    kmingen[i],
                    kmaxgen[i],
                    kminwav[i],
                    kmaxwav[i],
                    gen_intersection[i],
                    wav_intersection[i]);
// cleanup
            int kmingen_i, kmaxgen_i, kminwav_i, kmaxwav_i;
            bool gen_intersection_i, wav_intersection_i;
            get_onedim_intersections(intinfo[i],
                    lambda_j[i],
                    lambda_e[i],
                    lambda_k[i],
                    (((bc_[lambda_p][2*i])?0:2) + ((bc_[lambda_p][2*i+1])?0:1)),
                    min_type_i,
                    mu_j[i],
                    (((bc_[mu_p][2*i])?0:2) + ((bc_[mu_p][2*i+1])?0:1)),
                    kmingen_i,
                    kmaxgen_i,
                    kminwav_i,
                    kmaxwav_i,
                    gen_intersection_i,
                    wav_intersection_i);
// cleanup            
            assert (gen_intersection[i] == gen_intersection_i);
            assert (wav_intersection[i] == wav_intersection_i);
            if (gen_intersection_i)
            {
                assert (wav_intersection_i);
                assert (kmingen[i] == kmingen_i);
                assert (kmaxgen[i] == kmaxgen_i);
            }
            if (wav_intersection_i)
            {
                assert (kminwav[i] == kminwav_i);
                assert (kmaxwav[i] == kmaxwav_i);
            }
            
            if (!gen_intersection[i] && ! wav_intersection[i]) 
                return false;
        }
        return true;

    }


        
    template <class IBASIS, unsigned int DIM>
    Array1D<SampledMapping<DIM> > 
    QTBasis<IBASIS,DIM>::sampled_output(const unsigned int nunum,
            const double alpha,
            const bool primal,
            const int resolution)
    {
        Array1D<SampledMapping<DIM> > result;
        result.resize(get_nop());
        FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda (on nus patch)
        unsigned int nui_basisnum;
        Index nu(full_collection_[nunum]);
        MultiIndex<int,DIM> nu_j(nu.j()), nu_e(nu.e()), nu_k(nu.k());
        unsigned int nu_p (nu.p());
        
        for (unsigned int i = 0; i < DIM; i++)
        {
            nui_basisnum = (((bc_[nu_p][2*i])?0:2) + ((bc_[nu_p][2*i+1])?0:1));
            values[i] = evaluate(*bases_infact_[nui_basisnum],
                               typename IBASIS::Index(nu_j[i],
                                                      nu_e[i],
                                                      nu_k[i],
                                                      bases_infact_[nui_basisnum]),
                               primal,
                               resolution).values();
        }
        Point<DIM> a,b;
        for (unsigned int p=0; p < get_nop();++p)
        {
            for (unsigned int i=0; i<DIM; ++i)
            {
                a[i] = corners_[p][i];
                b[i] = a[i]+1;
            }
            Grid<DIM> grid(a,b, 1<<resolution);
            if (p == nu_p)
            {
                result[p] = SampledMapping<DIM>(a,b,values);
                result[p].mult(alpha);
            }
            else
            {
                result[p] = SampledMapping<DIM>(grid);
            }
        }
        MultiIndex<unsigned int, DIM> intinfo;
        //get_LMR_info(nunum, nu_j, nu_p, intinfo);
        get_LMR_info(nunum, nu_p, intinfo);
        unsigned int type;
        int centerpatchnumber;
        MultiIndex<bool,DIM> orientation;
         
        get_intersection_geometry(nu_p, intinfo,
                type, 
                centerpatchnumber, 
                orientation);
        Matrix<double> reflected_values(values[1].size(), values[0].size());
        Matrix<double> reflected2;
        int north, east, northeast;
        switch (type)
        {
            case 0:
                break;
            case 1:
                for (unsigned int m(0); m < reflected_values.row_dimension(); m++)
                    for (unsigned int n(0); n < reflected_values.column_dimension(); n++)
                        reflected_values(m,n) = values[1][m] * values[0][reflected_values.column_dimension() -1 -n];
                
                // copy values to the left or right neighbour
                // Observe: this works because SampledMapping.add ignores the grid!
                if (nu_p == centerpatchnumber)
                {
                    result[neighbours_[centerpatchnumber][1]].add(alpha, reflected_values);
                }
                else
                {
                    result[centerpatchnumber].add(alpha, reflected_values);
                }
                break;
            case 2:
                for (unsigned int m(0); m < reflected_values.row_dimension(); m++)
                    for (unsigned int n(0); n < reflected_values.column_dimension(); n++)
                        reflected_values(m,n) = values[1][reflected_values.row_dimension() -1 - m] * values[0][n];
                
                // copy values to the lower or upper neighbour
                // Observe: this works because SampledMapping.add ignores the grid!
                if (nu_p == centerpatchnumber)
                {
                    result[neighbours_[centerpatchnumber][3]].add(alpha, reflected_values);
                }
                else
                {
                    result[centerpatchnumber].add(alpha, reflected_values);
                }
                break;
            case 3:
                if (centerpatchnumber == -1)
                {
                    north = neighbours_[nu_p][0];
                    east = neighbours_[nu_p][2];
                }
                else
                {
                    north = neighbours_[centerpatchnumber][3]; 
                    east = neighbours_[centerpatchnumber][1];
                    if (north != -1)
                    {
                        northeast = neighbours_[north][1];
                    }
                    else
                    {
                        northeast = neighbours_[east][3];
                    }
                }
                // scaled version:
                for (unsigned int m(0); m < reflected_values.row_dimension(); m++)
                    for (unsigned int n(0); n < reflected_values.column_dimension(); n++)
                        reflected_values(m,n) = alpha * values[1][m] * values[0][reflected_values.column_dimension() -1 -n];
                // copy values to the neighbours
                // Observe: this works because SampledMapping.add ignores the grid!
                if (centerpatchnumber == nu_p)
                {
                    result[east].add(reflected_values );
                    reflected_values.reflect(reflected2);
                    result[north].add(reflected2) ;
                    if (northeast != -1)
                    {
                        // L-shaped?
                        result[nu_p].values().reflect(reflected2);
                        result[northeast].add(reflected2 );
                    }
                }
                else
                {
                    if (north == nu_p)
                    {
                        result[northeast].add(reflected_values);
                        reflected_values.reflect(reflected2);
                        result[centerpatchnumber].add(reflected2);
                        if (east != -1)
                        {
                            result[nu_p].values().reflect(reflected2);
                            result[east].add(reflected2);
                        }
                    }
                    else
                    {
                        if (east == nu_p)
                        {
                            result[centerpatchnumber].add(reflected_values);
                            reflected_values.reflect(reflected2);
                            result[northeast].add(reflected2);
                            if (north != -1)
                            {
                                result[nu_p].values().reflect(reflected2);
                                result[north].add(reflected2);
                            }
                        }
                        else
                        {
                            // northeast == nu_p
                            result[north].add(reflected_values);
                            reflected_values.reflect(reflected2);
                            result[east].add(reflected2);
                            if (centerpatchnumber != -1)
                            {
                                result[nu_p].values().reflect(reflected2);
                                result[centerpatchnumber].add(reflected2);
                            }
                        }
                    }
                }
                break;
            default:
                cout << "unrecognized geometry type! Maybe DIM==3??" << endl;
                abort();
                abort();
                break;
        }
        return result;
    }

    template <class IBASIS, unsigned int DIM>
    void 
    QTBasis<IBASIS,DIM>::sampled_output(Array1D<SampledMapping<DIM> > & previous_sampling,
            const unsigned int nunum,
            const double alpha,
            const bool primal,
            const int resolution)
    {
        FixedArray1D<Array1D<double>,DIM> values; // point values of the factors within psi_lambda (on nus patch)
        unsigned int nui_basisnum;
        Index nu(full_collection_[nunum]);
        MultiIndex<int,DIM> nu_j(nu.j()), nu_e(nu.e()), nu_k(nu.k());
        unsigned int nu_p (nu.p());
        
        for (unsigned int i = 0; i < DIM; i++)
        {
            nui_basisnum = (((bc_[nu_p][2*i])?0:2) + ((bc_[nu_p][2*i+1])?0:1));
            values[i] = evaluate(*bases_infact_[nui_basisnum],
                               typename IBASIS::Index(nu_j[i],
                                                      nu_e[i],
                                                      nu_k[i],
                                                      bases_infact_[nui_basisnum]),
                               primal,
                               resolution).values();
        }
        
        MultiIndex<unsigned int, DIM> intinfo;
        //get_LMR_info(nunum, nu_j, nu_p, intinfo);
        get_LMR_info(nunum, nu_p, intinfo);
        unsigned int type;
        int centerpatchnumber;
        MultiIndex<bool,DIM> orientation;
         
        get_intersection_geometry(nu_p, intinfo,
                type, 
                centerpatchnumber, 
                orientation);
        Matrix<double> matrix(values[1].size(), values[0].size());
        Matrix<double> matrix2;;
        for (unsigned int m(0); m < matrix.row_dimension(); m++)
            for (unsigned int n(0); n < matrix.column_dimension(); n++)
                matrix(m,n) = alpha * values[1][m] * values[0][n];
        previous_sampling[nu_p].add(matrix);
        int north, east,northeast;
        switch (type)
        {
            case 0:
                break;
            case 1:
                for (unsigned int m(0); m < matrix.row_dimension(); m++)
                    for (unsigned int n(0); n < matrix.column_dimension(); n++)
                        matrix(m,n) = values[1][m] * values[0][matrix.column_dimension() -1 -n];
                
                // copy values to the left or right neighbour
                // Observe: this works because SampledMapping.add ignores the grid!
                if (nu_p == centerpatchnumber)
                {
                    previous_sampling[neighbours_[centerpatchnumber][1]].add(alpha, matrix);
                }
                else
                {
                    previous_sampling[centerpatchnumber].add(alpha, matrix);
                }
                break;
            case 2:
                for (unsigned int m(0); m < matrix.row_dimension(); m++)
                    for (unsigned int n(0); n < matrix.column_dimension(); n++)
                        matrix(m,n) = values[1][matrix.row_dimension() -1 - m] * values[0][n];
                
                // copy values to the lower or upper neighbour
                // Observe: this works because SampledMapping.add ignores the grid!
                if (nu_p == centerpatchnumber)
                {
                    previous_sampling[neighbours_[centerpatchnumber][3]].add(alpha, matrix);
                }
                else
                {
                    previous_sampling[centerpatchnumber].add(alpha, matrix);
                }
                break;
            case 3:
                if (centerpatchnumber == -1)
                {
                    north = neighbours_[nu_p][0];
                    east = neighbours_[nu_p][2];
                }
                else
                {
                    north = neighbours_[centerpatchnumber][3];
                    east = neighbours_[centerpatchnumber][1];
                    if (north != -1)
                    {
                        northeast = neighbours_[north][1];
                    }
                    else
                    {
                        northeast = neighbours_[east][3];
                    }
                }
                if (centerpatchnumber == nu_p)
                {
                    if (northeast != -1)
                    {
                        matrix.reflect(matrix2);
                        previous_sampling[northeast].add(matrix2);
                    }
                    for (unsigned int m(0); m < matrix.row_dimension(); m++)
                        for (unsigned int n(0); n < matrix.column_dimension(); n++)
                            matrix(m,n) = alpha * values[1][m] * values[0][matrix.column_dimension() -1 -n];
                    previous_sampling[east].add(matrix);
                    matrix.reflect(matrix2);
                    previous_sampling[north].add(matrix2);
                }
                else
                {
                    if (north == nu_p)
                    {
                        if (east != -1)
                        {
                            matrix.reflect(matrix2);
                            previous_sampling[east].add(matrix2);
                        }
                        for (unsigned int m(0); m < matrix.row_dimension(); m++)
                            for (unsigned int n(0); n < matrix.column_dimension(); n++)
                                matrix(m,n) = alpha * values[1][m] * values[0][matrix.column_dimension() -1 -n];
                        previous_sampling[northeast].add(matrix);
                        matrix.reflect(matrix2);
                        previous_sampling[centerpatchnumber].add(matrix2);
                    }
                    else
                    {
                        if (east == nu_p)
                        {
                            if (north != -1)
                            {
                                matrix.reflect(matrix2);
                                previous_sampling[north].add(matrix2);
                            }
                            for (unsigned int m(0); m < matrix.row_dimension(); m++)
                                for (unsigned int n(0); n < matrix.column_dimension(); n++)
                                    matrix(m,n) = alpha * values[1][m] * values[0][matrix.column_dimension() -1 -n];
                            previous_sampling[centerpatchnumber].add(matrix);
                            matrix.reflect(matrix2);
                            previous_sampling[northeast].add(matrix2);
                        }
                        else
                        {
                            // northeast == nu_p
                            if (centerpatchnumber != -1)
                            {
                                matrix.reflect(matrix2);
                                previous_sampling[centerpatchnumber].add(matrix2);
                            }
                            for (unsigned int m(0); m < matrix.row_dimension(); m++)
                                for (unsigned int n(0); n < matrix.column_dimension(); n++)
                                    matrix(m,n) = alpha * values[1][m] * values[0][matrix.column_dimension() -1 -n];
                            previous_sampling[north].add(matrix);
                            matrix.reflect(matrix2);
                            previous_sampling[east].add(matrix2);
                        }
                    }
                }
                break;
            default:
                cout << "unrecognized geometry type!" << endl;
                cout << "method is only implemented for DIM = 2" << endl;
                abort();
                break;
        }
    }
    
    template <class IBASIS, unsigned int DIM>
    Array1D<SampledMapping<DIM> >
    QTBasis<IBASIS,DIM>::sampled_output(const InfiniteVector<double, int>& coeffs,
                const bool primal,
                const int resolution)
    {
        Array1D<SampledMapping<DIM> > result;
        InfiniteVector<double,int>::const_iterator it(coeffs.begin()), itend(coeffs.end());
        if (it != itend)
        {
            result = sampled_output(it.index(),*it,primal,resolution);
            ++it;
        }
        while (it != itend)
        {
            sampled_output(result, it.index(), *it, primal,resolution);
            ++it;
        }
        return result;
    }
    
    
        
}

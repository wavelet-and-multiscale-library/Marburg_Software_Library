// implementation for QTIndex.h

namespace WaveletTL
{
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>::QTIndex(const QTBASIS* basis)
    : basis_(basis), j_(), e_(), k_(), p_(0), num_(0)
    {}
    
#if 0
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>:: QTIndex(const level_type& j, const type_type& e, const translation_type& k, const int p, const QTBASIS* basis)
    : basis_(basis), j_(j), e_(e), k_(k), p_(p), num_(0)
    {
        abort();

        MathTL::FixedArray1D<std::map<int, int>,DIM> sizes; // store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
        int uptothislevel(0); // number of basis function below the current level
        int oncurrentlevel(1); // number of base elements on current level j
        int range(0); // determines whether a new entry has to be computed in "sizes"
        level_type currentlevel, j0(basis_->j0());
        type_type currenttype;
        for (unsigned int i=0; i < DIM; i++) {
            currentlevel[i]=j0[i];
            currenttype[i]=0;
            sizes[i][0] = basis_->bases()[i]->Deltasize(j0[i]); // Number of generators on level j0
            sizes[i][1] = basis_->bases()[i]->Nablasize(j0[i]); // Number of Wavelets on level j0
        }
        // iterate over all level up to j_ and add up number of basis functions up to the level
        if (currentlevel != j_ || currenttype != e_)
        while (true)
        {
            // compute number of functions on this level
            oncurrentlevel = 1;
            for (unsigned int i = 0; i < DIM; i++)
            {
                oncurrentlevel *= sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
            }
            uptothislevel += oncurrentlevel;
            // increase index = (currentlevel,currenttype)
            // "small loop"
            // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
            // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
            bool done = true;
            for (int i(DIM-1); i >= 0; i--)
            {
                // find first position on level j0
                if (currentlevel[i] == j0[i])
                {
                    if (currenttype[i] == 1)
                    {
                        currenttype[i]=0;
                    } else
                    {
                        currenttype[i]=1;
                        done = false;
                        break;
                    }
                }
            }
            // done == true means that all components with currentlevel[i]=j0[i] were wavelets.
            // the level j has to be increased.
            // iterate "big loop", meaning: "currentlevel++"

            if (done == true)
            {
            for (int i(DIM-1); i >= 0; i--)
            {
                if (i != 0)
                {
                    if (currentlevel[i] != j0[i])
                    {
                        // increase left neighbour
                        currentlevel[i-1]=currentlevel[i-1]+1;
                        if (currentlevel[i-1]-j0[i] == range) sizes[i-1][range+1]=basis ->bases()[i]->Nablasize(currentlevel[i-1]); // if needed compute and store new size information
                        currenttype[i-1]=1;
                        int temp = currentlevel[i]-j0[i];
                        currentlevel[i]=j0[i];
                        currenttype[i]=0;
                        currentlevel[DIM-1]=j0[DIM-1]+temp-1;
                        currenttype[DIM-1]= (temp == 1?0:1);
                        break;
                    }
                } else // i == 0. "big loop" arrived at the last index. We have to increase the level!
                {
                    range = range +1;
                    if (DIM == 1)
                    {
                        currenttype[i] = 1;
                        currentlevel[i]=currentlevel[i]+1;
                        sizes[0][range+1]=basis ->bases()[0]->Nablasize(currentlevel[0]); // if needed compute and store new size information
                    }
                    else
                    {
                        //currentlevel[DIM-1]=j0[DIM-1]+currentlevel[0]-j0[0]+1; currenttype[DIM-1]=1;
                        currentlevel[DIM-1]=j0[DIM-1]+range; currenttype[DIM-1]=1;
                        currenttype[0]=0; currentlevel[0]=j0[0];
                        sizes[DIM-1][range+1]=basis ->bases()[DIM-1]->Nablasize(currentlevel[DIM-1]); // if needed compute and store new size information
                    }
                break; // unnoetig, da i==0 gilt.
                }
            } // end of "big loop"
            }
            if (currentlevel == j_ && currenttype == e_) break;
        }
        // determine number of the function described by k in (currentlevel,currenttype)
        translation_type ktemp; // we only use that this is a multiindex
        for (unsigned int i = 0; i < DIM; i++)
            if (e_[i] == 0) ktemp[i]=k_[i]-basis_->bases()[i]->DeltaLmin();
            else ktemp[i]=k_[i]-basis_->bases()[i]->Nablamin();
        oncurrentlevel = ktemp[0];
        for (unsigned int i=1; i <= DIM-1; i++)
        {
            oncurrentlevel *= sizes[i][currentlevel[i]+currenttype[i]-j0[i]];
            oncurrentlevel += ktemp[i];
        }
        num_ =uptothislevel+oncurrentlevel;
    }
#endif
    
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>:: QTIndex(const level_type& j, const type_type& e, const translation_type& k, const int p, const int num, const QTBASIS* basis)
    : basis_(basis), j_(j), e_(e), k_(k), p_(p), num_(num)
    {
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS, DIM, QTBASIS>::QTIndex(const QTIndex& lambda)
    : basis_(lambda.basis_), j_(lambda.j_), e_(lambda.e_), k_(lambda.k_), p_(lambda.p_), num_(lambda.num_) {}

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS, DIM, QTBASIS>::QTIndex(const QTIndex* lambda)
    : basis_(lambda->basis_), j_(lambda->j_), e_(lambda->e_), k_(lambda->k_), p_(lambda->p_), num_(lambda->num_) {}

# if 0
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS, DIM, QTBASIS>::QTIndex(const int number, const QTBASIS* basis)
    : basis_(basis), num_(number)
    {
        abort();
        // very similar code to operator (++)
        // apply any changes here also there and vice versa
        Array1D<FixedArray1D<map<int, int>,DIM> > sizes; // Store number of basis elements. Generators on level j0 (0), wavelets on level j0 (1), j0+1 (2), ...
        sizes.resize(basis_->get_nop());
        int remains = number+1; // numbering begins at 0
        int oncurrentlevel(1); // number of base elements on current level j

        // int range(0); // number of level steps we have climbed so far. determines whether a new entry has to be computed in "sizes"

        //typename Index::type_type e;
        //typename Index::translation_type k;
        //level_type currentlevel;
        //type_type currenttype;
        //int currentpatch;
        Array1D<level_type> j0(basis_->j0());
        list<int> temp_list;
        int temp_int;

        for (unsigned int p=0; p < basis_->get_nop(); p++)
            {
                temp_list.push_back(p);
            }
        temp_list.sort(index_cmp<Array1D<level_type> > (j0,0));

        p_ = temp_list.front();
        temp_list.pop_front();

        j_ = j0[p_];

        for (unsigned int i=0; i < DIM; i++) 
        {
            e_[i]=0;
            for (unsigned int p=0; p< basis_->get_nop(); ++p)
            {
                sizes[p][i][0] = basis_->bases()[p][i]->Deltasize(j0[p][i]); // Number of generators on level j0
                sizes[p][i][1] = basis_->bases()[p][i]->Nablasize(j0[p][i]); // N o Wavelets
            }
            //oncurrentlevel *= sizes[p][i][0];
        }

        // iterate over all patches and levels. Add up number of basis functions till looked for level is reached
        //while (remains > oncurrentlevel)
        while(true)
        {
            // compute number of functions on this level
            oncurrentlevel = 1;
            for (unsigned int i = 0; i < DIM; i++)
            {
                oncurrentlevel *= sizes[p_][i][j_[i]+e_[i]-j0[p_][i]];
            }
            // CLEANUP
            if (remains < 10)
            {
                for (unsigned int i = 0; i < DIM; i++)
                {
                    // cout j,e,p,k :
                    cout << "i= " << i << " " << " sizes[p_][i][j_[i]+e_[i]-j0[p_][i]] = " << sizes[p_][i][j_[i]+e_[i]-j0[p_][i]] << endl;
                }
            }
            if (remains > oncurrentlevel)
            {
                //substract number of basis functions on current_index=(j_,e_,p_) and increase current_index
                remains -= oncurrentlevel;
            }
            else
            {
                // break if we are at the right level
                break;
            }

            // "small loop" "e_++"
            // iterate over all combinations of generators/wavelets in all dimensions on the current level
            // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)

            bool done = false;
            for (int i(DIM-1); i >= 0; i--)
            {
                // find first position on level j0
                if (j_[i] == j0[p_][i])
                {
                    if (e_[i] == 1)
                    {
                        e_[i]=0;
                    } else
                    {
                        e_[i]=1;
                        done = true; // ++e_ was successful
                        break;
                    }
                }
            }
            if (done)
            {
                continue;
            }
            // done == false means that all components with j_[i]=j0[p][i] were wavelets.
            // "p_++"
            while ((p_ < basis_->get_nop()-1) && !done)
            {
                p_++;
                done = true;
                temp_int = 0;
                // check if level j_ exists for the current patch.
                // that is not the same as using > or lex from MultiIndex
                while (done == true)
                {
                    done = (j_[temp_int] >= j0[p_][temp_int]);
                    temp_int++;
                    if (temp_int == DIM)
                    {
                        break;
                    }
                }
                if (done == true) // level is valid for p_
                {
                    for (unsigned int i=0; i<DIM; i++)
                    {
                        temp_int = j_[i]-j0[p_][i];
                        // minimal level may have changed
                        if (temp_int == 0)
                        {
                            e_[i]= 0;
                        }
                        else
                        {
                            e_[i] = 1;
                            if (sizes[p_][i].size()<temp_int+2) // if needed compute and store new size information
                            {
                                // compute sizes:
                                sizes[p_][i][temp_int+2] = basis_->bases()[p_][i]->Nablasize(j_[i]);
                            }
                        }
                    }
                    //return *this;
                }
            }

            if (done)
            {
                continue;
            }

            // If it is impossible to increase p_ for the current level j_,
            // then we have to increase the norm of j_ by 1.


            // Special case 1D:
            // it is considered s.t. we can asert DIM > 1 for the rest of the method

            if (DIM == 1)
            {
                j_[0]= j_[0]+1;
                // Asumption: the minimal levels on different patches differ at most by 1
                // Thus the first patch allows the level j_
                p_ = 0;
                e_[0] = (j_[0] != j0[0][0]);

                temp_int = j_[0]-j0[0][0];
                if (sizes[0][0].size()<temp_int+2) // if needed compute and store new size information
                {
                    // compute sizes:
                    sizes[0][0][temp_int+2] = basis_->bases()[0][0]->Nablasize(j_[0]);
                }
                // return *this;

                // Use this code, to allow arbitrary minimal levels:
                /*
                for (unsigned int p=0; p < basis_->get_nop(); ++p)
                {
                    if (j_[0] >= j0[p][0])
                    {
                        p_ = p;
                        e_[0] = (j_[0] != j0[p][0]);
                 * // compute sizes:
                        //return *this;
                    }
                }
                */
                continue;
            }

            // try to increase the sublevel,
            // that is increasing j_ without increasing its 1-norm.

            FixedArray1D<list<int>*, DIM-1> incrementable; // store the patches which allow for an increment
            Array1D<int> min_suffix_norms; // for each patch store the norm of the suffix (from current to last index)
            min_suffix_norms.resize(basis_->get_nop());
            temp_list.clear(); // store the patches which allow the prefix up to the current dimension

            incrementable[0]=new list<int>;

            for (unsigned int p=0; p < basis_->get_nop(); p++)
            {
                if ( (j0[p][0] <= j_[0]+1) && (j0[p][1] < j_[1]) ) // check if decrement at pos=1 and increment at pos =0 are possible
                {
                    incrementable[0]->push_back(p);
                }
                // CLEANUP
                if (remains < 10)
                {
                    cout << "j0[p][0] = " << j0[p][0] << " j_[0] = " << j_[0] << " j0[p][1] = " << j0[p][1] << " j_[1] = " << j_[1] << endl;
                }
            }
            
            if (DIM > 2)
            {
                for (unsigned int p=0; p < basis_->get_nop(); p++)
                {
                    if (j0[p][0] <= j_[0])
                    {
                        temp_list.push_back(p);
                    }
                }
            }

            for (unsigned int i=1; i< DIM-1; i++)
            {
                incrementable[i]=new list<int>;
                for (list<int>::iterator it(temp_list.begin()), itend(temp_list.end()); it != itend;)
                {
                    if ( (j0[*it][i] <= j_[i]+1) && (j0[*it][i+1] < j_[i+1]))
                    {
                        incrementable[i]->push_back(*it);
                    }
                    if (j0[*it][i] > j_[i])
                    {
                        temp_list.erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }
            }
            // from right to left: search for an index i to decrease j_[i]

            // i == DIM-2
            if (!incrementable[DIM-2]->empty())
            {
                // increment at DIM-2, decrement at DIM-1 => no suffix.
                // sorting increment[DIM-2] by patchnumber is sufficient
                incrementable[DIM-2]->sort();
                //accept!
                p_ = incrementable[DIM-2]->front();
                j_[DIM-2] += 1;
                e_[DIM-2] = (j_[DIM-2] != j0[p_][DIM-2]);
                j_[DIM-1] -= 1;
                e_[DIM-1] = (j_[DIM-1] != j0[p_][DIM-1]);
                temp_int = j_[DIM-2]-j0[p_][DIM-2];
                if (sizes[p_][DIM-2].size()<temp_int+2) // if needed compute and store new size information
                {
                    sizes[p_][DIM-2][temp_int+2] = basis_->bases()[p_][DIM-2]->Nablasize(j_[DIM-2]);
                }

                temp_int = j_[DIM-1]-j0[p_][DIM-1];
                if (sizes[p_][DIM-1].size()<temp_int+2) // if needed compute and store new size information
                {
                    // I think this line is only executed if minimal levels differ by more than one!
                    sizes[p_][DIM-1][temp_int+2] = basis_->bases()[p_][DIM-1]->Nablasize(j_[DIM-1]);
                }
                continue;
                //return *this;
            }
            else
            {
                temp_int = j_[DIM-1]; // store j_[i+1] + ... + j_[DIM-1]
                for (unsigned int p = 0; p < basis_->get_nop();++p)
                {
                    min_suffix_norms[p]=basis_->j0()[p][DIM-1];
                }
            }

            int temp_int2;
            for (int i=DIM-3; i>=0; --i)
            {
                temp_int += j_[i+1];
                for (unsigned int p = 0; p < basis_->get_nop();++p)
                {
                    min_suffix_norms[p]+=j0[p][i+1];
                }
                if (!incrementable[i]->empty())
                {
                    // sort increment set according to j0[p][i+1,...,DIM-1]
                    incrementable[i]->sort(index_cmp_ignoreLastEntry<Array1D<MultiIndex<int,DIM> > > (j0,i+1));

                    for (list<int>::const_iterator it(incrementable[i]->begin()), itend(incrementable[i]->end()); it != itend; ++it)
                    {
                        // accept if j_[i+1] + ... + j_[DIM-1] > j0[p][i+1] + ... j0[p][DIM-1]
                        if (min_suffix_norms[*it] < temp_int)
                        {
                            //accept!
                            p_ = *it;
                            j_[i] = j_[i]+1;
                            e_[i] = (j_[i] != j0[p_][i]);
                            for ( unsigned int j = i+1; j <= DIM -2; ++j)
                            {
                                j_[j] = j0[p_][j];
                                e_[j] = 0;
                            }
                            j_[DIM-1] = j0[p_][DIM-1] + temp_int - min_suffix_norms[p_] -1;
                            e_[DIM-1] = (j_[DIM-1] != j0[p_][DIM-1]);

                            temp_int2 = j_[i]-j0[p_][i];
                            if (sizes[p_][i].size()<temp_int2+2) // if needed compute and store new size information
                            {
                                sizes[p_][i][temp_int2+2] = basis_->bases()[p_][i]->Nablasize(j_[i]);
                            }
                            //temp_int2 = temp_int - min_suffix_norms[p_] -1;
                            if (sizes[p_][DIM-1].size()<temp_int - min_suffix_norms[p_] +1) // if needed compute and store new size information
                            {
                                // I think this line is only executed if minimal levels differ by more than one!
                                sizes[p_][DIM-1][temp_int - min_suffix_norms[p_] +1] = basis_->bases()[p_][DIM-1]->Nablasize(j_[DIM-1]);
                            }
                            done = true;
                            break;
                            //return *this;
                        }
                    }
                }
                if (done == true)
                {
                    break;
                }
            }

            if (done == true)
            {
                continue;
            }

            // "big loop" : increase j_ such that the 1-norm increases by 1
            temp_int = multi_degree(j_)+1; // norm of the new level
            temp_list.clear();
            for (unsigned int p=0; p< basis_->get_nop(); ++p)
            {
                temp_list.push_back(p);
            }
            temp_list.sort(index_cmp_ignoreLastEntry< Array1D<MultiIndex<int,DIM> > > (j0,0));

            for (list<int>::const_iterator it(temp_list.begin()),itend(temp_list.end());it!=itend;++it)
            {
                temp_int2 = multi_degree(j0[*it]);
                if ( temp_int2 <= temp_int)
                {
                    //accept this patch
                    j_ = j0[*it];
                    p_= (*it);
                    for (unsigned int i = 0; i<DIM-1; ++i)
                    {
                        e_[i]=0;
                    }
                    j_[DIM-1] = j0[p_][DIM-1] + temp_int - temp_int2;
                    e_[DIM-1] = (j_[DIM-1] != j0[p_][DIM-1]);

                    sizes[p_][DIM-1][temp_int - temp_int2] = basis_->bases()[p_][DIM-1]->Nablasize(j_[DIM-1]);
                    // I think the computation of size here is always needed if minlevels differ less than 1
                    assert (temp_int - temp_int2 > 0);
                    // done = true;
                    break;
                    //return *this;
                }
            }
        } // end of while

        // determine k corresponding to the number of the basis function given by "remains" (on level (j_,e_) )
        unsigned int modul;
        //j_ = currentlevel;
        //e_ = currenttype;
        remains -= 1; // numbering begins at 0
        for (int i = DIM-1; i > 0; i--)
        {
                modul = sizes[p_][i][j_[i]+e_[i]-j0[p_][i]];
                k_[i]= remains % modul+(e_[i] == 0 ? basis->bases()[p_][i]->DeltaLmin():basis->bases()[p_][i]->Nablamin());
                remains = remains / modul;
        }
        k_[0]=remains +(e_[0] == 0 ? basis->bases()[p_][0]->DeltaLmin():basis->bases()[p_][0]->Nablamin());
    }
#endif

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>&
    QTIndex<IBASIS,DIM,QTBASIS>::operator = (const QTIndex<IBASIS,DIM,QTBASIS>& lambda)
    {
    j_ = lambda.j();
    e_ = lambda.e();
    k_ = lambda.k();
    p_ = lambda.p();
    basis_ = lambda.basis();
    num_ = lambda.number();
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    bool
    QTIndex<IBASIS, DIM, QTBASIS>::operator == (const QTIndex<IBASIS, DIM, QTBASIS>& lambda) const
    {
        return (p_ == lambda.p() &&
                j_ == lambda.j() &&
                e_ == lambda.e() &&
                k_ == lambda.k() &&
                num_ == lambda.number());
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>&
    QTIndex<IBASIS,DIM,QTBASIS>::operator ++ ()
    {
        if (num_ == -1) return *this;
        Array1D<level_type> j0(basis_->j0());
        num_++;
        // try to increase translation index k
        bool done = true;
        int temp_int;
        for (int i = DIM-1; i >= 0; i--)
        {
            // store number of last index in current dimension
            temp_int = (e_[i] == 0 ? basis_->bases()[p_][i]->DeltaRmax(j_[i])
                                   : basis_->bases()[p_][i]->Nablamax(j_[i]));
            if (k_[i] == temp_int)
            {
                k_[i] = (e_[i] == 0 ? basis_->bases()[p_][i]->DeltaLmin()
                                    : basis_->bases()[p_][i]->Nablamin());
                done = (i != 0);
                //jplusplus = (i == 0);
            } else
            {
                ++k_[i];
                break;
            }
        }

        if (done == true) return *this;

        // try to increase type e
        // (k_[i] is set to the minimal value DeltaLmin or Nablamin)
        // iterate over all combinations of generators/wavelets for all dimensions with j_[p][i]=j0[p][i]
        // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)

        for (int i(DIM-1); i >= 0; i--)
        {
            // find first position on level j0
            if (j_[i] == j0[p_][i])
            {
                if (e_[i] == 1)
                {
                    e_[i]=0;
                    k_[i]=basis_->bases()[p_][i]->DeltaLmin();
                } else
                {
                    e_[i]=1;
                    k_[i]=basis_->bases()[p_][i]->Nablamin();
                    done = true;
                    break;
                }
            }
        }

        if (done == true) return *this;

        // try to increase patch p
        while (p_ < basis_->get_nop()-1)
        {
            p_++;
            done = true;
            temp_int = 0;
            // check if level j_ exists for the current patch.
            // that is not the same as using > or lex from MultiIndex
            while (done == true)
            {
                done = (j_[temp_int] >= j0[p_][temp_int]);
                temp_int++;
                if (temp_int == DIM)
                {
                    break;
                }
            }
            if (done == true)
            {
                for (unsigned int i=0; i<DIM; i++)
                {
                    // minimal level may have changed
                    e_[i] = (j_[i] != j0[p_][i]);
                    k_[i] = (e_[i] == 0 ? basis_->bases()[p_][i]->DeltaLmin()
                                        : basis_->bases()[p_][i]->Nablamin());
                }
                return *this;
            }
        }

        // Special case 1D:
        // If it is impossible to increase p_ for the current level j_,
        // then we have to increase the norm of j_ by 1. This case is considered
        // here s.t. we can assert DIM > 1 for the rest of the method
        if (DIM == 1)
        {
            j_[0]= j_[0]+1;
            // Assumption: the minimal levels on different patches differ at most by 1
            p_ = 0;
            e_[0] = (j_[0] != j0[0][0]);
            k_[0] = (e_[0] == 0 ? basis_->bases()[0][0]->DeltaLmin()
                                : basis_->bases()[0][0]->Nablamin());
            return *this;
            // Use this code, to allow arbitrary minimal levels:
            /*
            for (unsigned int p=0; p < basis_->get_nop(); ++p)
            {
                if (j_[0] >= j0[p][0])
                {
                    p_ = p;
                    e_[0] = (j_[0] != j0[p_][0]);
                    k_[0] = (e_[0] == 0 ? basis_->bases()[p_][0]->DeltaLmin()
                                        : basis_->bases()[p_][0]->Nablamin());
                    return *this;
                }
            }
            */
        }
        assert (DIM > 1);

        // try to increase the sublevel,
        // that is increasing j_ without increasing its 1-norm.
        // A strategy that allows for arbitrary minimal levels would be:
        // Find the first level on each patch that comes after j_ and take of all lowest ones the one with the least patchnumber.
        // Not going to do this. Assume: The leftmost index that has to be increased is only increased by one.
        // That is the same as assuming that the minimal levels on each patch and in each direction differ at most by 1.
        // 
        // Example: minimal levels 223 and 422. j_=223. The third position in j_[2] cannot be
        // reduced even though there is a patch with j0[2] < j_[2]. We need to check which
        // patches allow the prefix of the third entry "23" ="2 (2+1) " and search this set for a patch
        // which allows the decrement of the current position.
        // Example: minlevels 223, 214. j_=223. j_[2] cannot be reduced. j_[1] can,
        // but the resulting level 313 is illegal. This means that we not only need to check
        // the prefix and if the current positions allows for an decrement, but also if the
        // resulting suffix is legal. There may be different possible suffices, we select the lowest possible.
        // Example: minlevels 222, 232. j_=252. j_[2] cannot be reduced. j_[1] can and (j_[0]+1) is legal.
        // We need to give priority to the patch which allows for the biggest reduction!



        FixedArray1D<std::list<int>*, DIM-1> incrementable; // store the patches which allow for an increment
        Array1D<int> min_suffix_norms; // for each patch store the norm of the suffix (from current to last index)
        min_suffix_norms.resize(basis_->get_nop());
        //int min_suffix_norm;
        
        //FixedArray1D<list<int>*, DIM-1> prefix;
        list<int>* temp_list;

        incrementable[0]=new list<int>;
        temp_list = new list<int>; // store the patches which allow the prefix up to the current dimension

        for (unsigned int p=0; p < basis_->get_nop(); p++)
        {
            if ( (j0[p][0] <= j_[0]+1) && (j0[p][1] < j_[1]) ) // check if decrement at pos=1 and increment at pos =0 are possible
            {
                incrementable[0]->push_back(p);
            }
        }
        if (DIM > 2)
        {
            for (unsigned int p=0; p < basis_->get_nop(); p++)
            {
                if (j0[p][0] <= j_[0])
                {
                    temp_list->push_back(p);
                }
            }
            //incrementable[0]->sort(index_cmp< map <int,MultiIndex<int,DIM>&> (basis_->j0(),2)); // ersetze map<...> durch den korrekten Rückgabewert
            //incrementable[0]->sort(index_cmp< Array1D<MultiIndex<int,DIM> >&> (basis_->j0(),2) ); // ersetze map<...> durch den korrekten Rückgabewert

        }

        for (unsigned int i=1; i< DIM-1; i++)
        {
            incrementable[i]=new list<int>;
            for (list<int>::iterator it(temp_list->begin()), itend(temp_list->end()); it != itend;)
            {
                if ( (j0[*it][i] <= j_[i]+1) && (j0[*it][i+1] < j_[i+1]))
                {
                    incrementable[i]->push_back(*it);
                    //if (DIM > 2+i)
                    //{
                    //    // no sorting needed in the last iteration. Ordering by patchnumber is sufficient in this case.
                    //    //incrementable[i]->sort(index_cmp<map<int,MultiIndex<int,DIM>&>(basis_->j0(),2+i)); // ersetze map<...> durch den korrekten Rückgabewert
                    //    incrementable[i]->sort(index_cmp< Array1D<MultiIndex<int,DIM> >&> (basis_->j0(),2+i) ); // ersetze map<...> durch den korrekten Rückgabewert
                    //}
                }
                if (j0[*it][i] > j_[i])
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
            // increment at DIM-2, decrement at DIM-1 => no suffix.
            // sorting increment[DIM-2] by patchnumber is sufficient
            incrementable[DIM-2]->sort();
            //accept!
            // only works if minimal differ only by 1
            p_ = incrementable[DIM-2]->front();
            j_[DIM-2] += 1;
            e_[DIM-2] = (j_[DIM-2] != j0[p_][DIM-2]);
            k_[DIM-2] = (e_[DIM-2] == 0 ? basis_->bases()[p_][DIM-2]->DeltaLmin()
                                        : basis_->bases()[p_][DIM-2]->Nablamin());
            j_[DIM-1] -= 1;
            e_[DIM-1] = (j_[DIM-1] != j0[p_][DIM-1]);
            k_[DIM-1] = (e_[DIM-1] == 0 ? basis_->bases()[p_][DIM-1]->DeltaLmin()
                                        : basis_->bases()[p_][DIM-1]->Nablamin());
            return *this;
        }
        temp_int = j_[DIM-1]; // store j_[i+1] + ... + j_[DIM-1]
        for (unsigned int p = 0; p < basis_->get_nop();++p)
        {
            min_suffix_norms[p]=basis_->j0()[p][DIM-1];
        }

        for (int i=DIM-3; i>=0; --i)
        {
            temp_int += j_[i+1];
            // The "incrementable" sets are not nested in general!
            // 2 possibilities: compute min_suffix_norms for all patches and lengths,
            // or compute the value as it is needed. The former stores more values but 
            // sometimes we may reduce the "incrementable" set. In the latter case we 
            // always need to sort the whole set, but only for the first few values 
            // minsuffix values need to be computed.
            for (unsigned int p = 0; p < basis_->get_nop();++p)
            {
                min_suffix_norms[p]+=basis_->j0()[p][i+1];
            }
            if (!incrementable[i]->empty())
            {
                // we want to increment at position i.
                // the new level j_ will satisfy j_[k] = j0()[p][k] for k = i+1,...,DIM-2
                // and j_[DIM-1] >= j0()[p][DIM-1]
                // This means: we may deduce the patch p by sorting the set incrementable[i] w.r.t. 
                // j0()[i+1,...,DIM-2] and then by patchnumber. 
                // Observe that j0()[?][DIM-1] is not relevant for the sorting. 
                // However, this value is relevant in the sense, that the new level j_ needs to be legal. We may not increase \|j_\| by accident.
                // Example: j0() = [4222, 4223, 3322]. j_ = 3322. 
                // Incrementable[2] and [1] are empty. incrementable[0] contains patch 0,1.
                // Increment of j_ at pos 0, decrement at pos 1 yields 4222. 
                // Observe that patch 1 has j0()[1][DIM-1] == 3. Increment of j_ at pos 0, decrement at pos 1 and setting the rest to j0()[1][i+1,..,DIM-1] would yield 
                // 4223. This would increase the norm of j_ which we do not want! 
                // We therefore need to check in incrementable, whether the norm of a suffix j0()[p][i+1,...,DIM-1] is small enough!
                
                // sort increment set according to j0[p][i+1,...,DIM-1] (wrong!)
                //incrementable[i]->sort(index_cmp<Array1D<MultiIndex<int,DIM> > > (basis_->j0(),i+1));
                
                // sort increment set according to j0[p][i+1,...,DIM-2]
                incrementable[i]->sort(index_cmp_ignoreLastEntry<Array1D<MultiIndex<int,DIM> > > (basis_->j0(),i+1));
                
                
                //check min_suffix_norm condition
                for (list<int>::const_iterator it(incrementable[i]->begin()), itend(incrementable[i]->end()); it != itend; ++it)
                {
                    // accept if j_[i+1] + ... + j_[DIM-1] > j0[p][i+1] + ... j0[p][DIM-1]

                    /* // compute min_suffix_norms as they are needed:
                    min_suffix_norm = 0;
                    for (unsigned int j= i+1; j <= DIM-1; j++)
                    {
                        //min_suffix_norms[*it] = suffix_norms[*it] + j0[*it][j];
                        min_suffix_norm += j0[*it][j];
                    }
                     */
                    if (min_suffix_norms[*it] < temp_int)
                    {
                        //accept!
                        // works for arbitrary minimal levels
                        p_ = *it;
                        j_[i] = j_[i]+1;
                        e_[i] = (j_[i] != j0[p_][i]);
                        k_[i] = (e_[i] == 0 ? basis_->bases()[p_][i]->DeltaLmin()
                                            : basis_->bases()[p_][i]->Nablamin());
                        for ( unsigned int j = i+1; j <= DIM -2; ++j)
                        {
                            j_[j] = j0[p_][j];
                            e_[j] = 0;
                            k_[j] = basis_->bases()[p_][j]->DeltaLmin();
                        }
                        j_[DIM-1] = j0[p_][DIM-1] + temp_int - min_suffix_norms[p_] -1;
                        e_[DIM-1] = (j_[DIM-1] != j0[p_][DIM-1]);
                        k_[DIM-1] = (e_[DIM-1] == 0 ? basis_->bases()[p_][DIM-1]->DeltaLmin()
                                                    : basis_->bases()[p_][DIM-1]->Nablamin());
                        return *this;
                    }
                }
            }
        }

        // "big loop" : increase j_ such that the 1-norm increases by 1

        // hereafter stands similar code as in first_generator()/wavelet(unsigned int j)
        assert ((temp_int + j_[0]) == multi_degree(j_));
        for (unsigned int p = 0; p < basis_->get_nop();++p)
        {
            min_suffix_norms[p]+=j0[p][0];
            assert (min_suffix_norms[p] == multi_degree(j0[p]));
        }
        temp_int += j_[0]+1; // norm of the new level
        //temp_int = multi_degree(j_)+1; // norm of the new level
        temp_list->clear();
        for (unsigned int p=0; p< basis_->get_nop(); p++)
        {
            temp_list->push_back(p);
        }
        // We deduce the new patchnumber by sorting the minimal levels
        temp_list-> sort(index_cmp_ignoreLastEntry< Array1D<MultiIndex<int,DIM> > > (j0,0));
        // The following code is checks again the min_suffix_condition as above. 
        // This is only relevant if the minimal levels of patches are allowed to differ by MORE than 1. 
        // If they do not, then the code should always select the first entry
        for (list<int>::const_iterator it(temp_list->begin()),itend(temp_list->end());it!=itend;++it)
        {
            if (min_suffix_norms[*it] <= temp_int)
            {
                //accept this patch
                j_ = j0[*it];
                p_= (*it);
                for (unsigned int i = 0; i<DIM-1; ++i)
                {
                    e_[i]=0;
                    k_[i] = basis_->bases()[p_][i]->DeltaLmin();
                }
                j_[DIM-1] = j0[p_][DIM-1] + temp_int - min_suffix_norms[*it];
                e_[DIM-1] = (j_[DIM-1] != j0[p_][DIM-1]);
                k_[DIM-1] = (e_[DIM-1] == 0 ? basis_->bases()[p_][DIM-1]->DeltaLmin()
                                            : basis_->bases()[p_][DIM-1]->Nablamin());
                return *this;
            }
            else
            {
                // since we assumed the minimal levels differ by at most 1
                abort();
            }
        }
        abort();
        return *this;
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    bool
    QTIndex<IBASIS,DIM,QTBASIS>::operator < (const QTIndex& lambda) const
    {
        // Ordering primary by level j as in MultiIndex
        // (ordering of \N^dim, that is the distance from 0, that is the same as
        // ordering first by 1-norm of j and in the case of equal norms lexicographical in j.)
        // then lexicographical in p,e,k.
        return (multi_degree(j_) < multi_degree(lambda.j()) ||
                (multi_degree(j_) == multi_degree(lambda.j()) &&
//                     (j_ < lambda.j() ||
                 (j_.lex(lambda.j()) || // should be a tiny bit faster
                  (j_ == lambda.j() &&
                   (p_ < lambda.p() ||
                    (p_ == lambda.p() &&
                     (e_.lex(lambda.e()) ||
                      (e_ == lambda.e() && k_.lex(lambda.k()))
                     )
                    )
                   )
                  )
                 )
                )
               );
    }

#if 0 //first / last routines
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>
    first_generator(const QTBASIS* basis)
    {
        return QTIndex<IBASIS,DIM,QTBASIS>(0, basis);
    }


    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>
    last_generator(const QTBASIS* basis)
    {
            typename QTIndex<IBASIS,DIM,QTBASIS>::level_type j;
            typename QTIndex<IBASIS,DIM,QTBASIS>::type_type e;
            typename QTIndex<IBASIS,DIM,QTBASIS>::translation_type k;
            j = basis->j0();
            for (unsigned int i = 0; i < DIM; i++)
            {
                    e[i] = 0;
                    k[i] = basis->bases()[i]->DeltaRmax(j[i]);
            }
            return QTIndex<IBASIS,DIM,QTBASIS>(j, e, k, basis);
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>
    first_wavelet(const QTBASIS* basis, const typename QTIndex<IBASIS,DIM,QTBASIS>::level_type j)
    {
        assert(multi_degree(j) > multi_degree(basis->j0())
                || ((multi_degree(j) == multi_degree(basis->j0())) && (basis->j0() <= j))
              );
        typename QTIndex<IBASIS,DIM,QTBASIS>::type_type e;
        typename QTIndex<IBASIS,DIM,QTBASIS>::translation_type k;
        bool first_level = true;
        for (unsigned int i = 0; i < DIM; i++) {
            if (j[i] == basis->bases()[i]->j0())
            {
                e[i] = 0;
                k[i] = basis->bases()[i]->DeltaLmin();
            } else
            {
                e[i] = 1;
                k[i] = basis->bases()[i]->Nablamin();
                first_level = false;
            }
        }
        if (first_level == true)
        {
            e[DIM-1] = 1;
            k[DIM-1] = basis->bases()[DIM-1]->Nablamin();
        }
        return QTIndex<IBASIS,DIM,QTBASIS>(j, e, k, basis);
    }


    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>
    first_wavelet(const QTBASIS* basis, const int level)
    {
        assert(level >= multi_degree(basis->j0()) );
        typename QTIndex<IBASIS,DIM,QTBASIS>::level_type j;
        typename QTIndex<IBASIS,DIM,QTBASIS>::type_type e;
        typename QTIndex<IBASIS,DIM,QTBASIS>::translation_type k;

        j[DIM-1] = basis->jo()[DIM-1] + level - multi_degree(basis->j0());
        e[DIM-1] = 1;
        k[DIM-1] = basis->basis()[DIM-1]->Nablamin();

        for (unsigned int i = 0; i < DIM-1; i++)
        {
            j[i] = basis->j0()[i];
            e[i] = 0;
            k[i] = basis->bases[i]->DeltaLmin();
        }
        return QTIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, basis);
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>
    last_wavelet(const QTBASIS* basis, const typename QTIndex<IBASIS,DIM,QTBASIS>::level_type j)
    {
        assert(multi_degree(j) > multi_degree(basis->j0())
                || ((multi_degree(j) == multi_degree(basis->j0())) && (basis->j0() <= j))
              );
        typename QTIndex<IBASIS,DIM,QTBASIS>::type_type e;
        typename QTIndex<IBASIS,DIM,QTBASIS>::translation_type k;
        for (unsigned int i = 0; i < DIM; i++) {
            e[i] = 1;
            k[i] = basis->bases()[i]->Nablamax(j[i]);
        }
        return QTIndex<IBASIS,DIM,QTBASIS>(j, e, k, basis);
    }

    // PERFORMANCE !! viele Aufrufe von j0()
    template <class IBASIS, unsigned int DIM, class QTBASIS>
    QTIndex<IBASIS,DIM,QTBASIS>
    last_wavelet(const QTBASIS* basis, const unsigned int level)
    {
        assert(level >= multi_degree(basis->j0()) );
        typename QTIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::level_type j;
        typename QTIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::type_type e;
        typename QTIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >::translation_type k;
        j[0]= basis->j0()[0]+level-multi_degree(basis->j0());
        e[0]=1;
        k[0]= basis->bases()[0]->Nablamax(j[0]);
        for (unsigned int i = 1; i < DIM; i++) {
            j[i] = basis->j0()[i];
            e[i] = 1;
            k[i] = basis->bases()[i]->Nablamax(basis->bases()[i]->j0());
        }
        return QTIndex<IBASIS,DIM,TensorBasis<IBASIS,DIM> >(j, e, k, basis);
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    const int
    first_generator_num(const QTBASIS* basis)
    {
            return 0;
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    const int
    last_generator_num(const QTBASIS* basis)
    {
            int res=1;
            for (unsigned int i = 0; i < DIM; i++)
                    res *= basis->bases()[i]->Deltasize(basis->bases()[i]->j0());
            return res-1;
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    const int
    first_wavelet_num(const QTBASIS* basis, const typename QTIndex<IBASIS,DIM,QTBASIS>::level_type j)
    {
        // Assertions are checked in first_wavelet
        // for (unsigned int i = 0; i < DIM; i++)
        //     assert(j[i] >= (basis->bases()[i]->j0()));
        QTIndex<IBASIS,DIM,QTBASIS> temp (first_wavelet<IBASIS,DIM,QTBASIS>(basis,j));
        return temp.number();
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    const int
    last_wavelet_num(const QTBASIS* basis, const typename QTIndex<IBASIS,DIM,QTBASIS>::level_type j)
    {
        // Assertions are checked in last_wavelet
        // for (unsigned int i = 0; i < DIM; i++)
        //     assert(j[i] >= (basis->bases()[i]->j0()));
        QTIndex<IBASIS,DIM,QTBASIS> temp (last_wavelet<IBASIS,DIM,QTBASIS>(basis,j));
        return temp.number();
    }

    template <class IBASIS, unsigned int DIM, class QTBASIS>
    const int
    last_wavelet_num(const QTBASIS* basis, const unsigned int level)
    {
        // Assertions are checked in last_wavelet
        // for (unsigned int i = 0; i < DIM; i++)
        //     assert(j[i] >= (basis->bases()[i]->j0()));
        QTIndex<IBASIS,DIM,QTBASIS> temp (last_wavelet<IBASIS,DIM,QTBASIS>(basis,level));
        return temp.number();
    }
#endif
}


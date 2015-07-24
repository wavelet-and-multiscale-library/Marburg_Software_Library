// implementation for tbasis_support.h

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>
#include <iostream>

using MathTL::multi_degree;
using MathTL::FixedArray1D;

namespace WaveletTL
{
    template <class IBASIS, unsigned int DIM>
    inline
    void
    support(const TensorBasis<IBASIS,DIM>& basis,
            const typename TensorBasis<IBASIS,DIM>::Index& lambda,
            typename TensorBasis<IBASIS,DIM>::Support& supp)
    {
        basis.support(lambda, supp);
    }

    template <class IBASIS, unsigned int DIM>
    bool
    intersect_supports(const TensorBasis<IBASIS,DIM>& basis,
                       const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                       const typename TensorBasis<IBASIS,DIM>::Index& mu,
                       typename TensorBasis<IBASIS,DIM>::Support& supp)
    {
        typename TensorBasis<IBASIS,DIM>::Support supp_lambda;
        basis.support(lambda, supp_lambda);
        typename TensorBasis<IBASIS,DIM>::Support supp_mu;
        basis.support(mu, supp_mu);
        // determine support intersection granularity,
        // adjust single support granularities if necessary
        for (unsigned int i=0;i<DIM;i++)
        {
            if (supp_lambda.j[i] > supp_mu.j[i]) {
                supp.j[i] = supp_lambda.j[i];
                const int adjust = 1<<(supp_lambda.j[i]-supp_mu.j[i]);
                supp_mu.a[i] *= adjust;
                supp_mu.b[i] *= adjust;
            } else {
                supp.j[i] = supp_mu.j[i];
                const int adjust = 1<<(supp_mu.j[i]-supp_lambda.j[i]);
                supp_lambda.a[i] *= adjust;
                supp_lambda.b[i] *= adjust;
            }
            supp.a[i] = std::max(supp_lambda.a[i],supp_mu.a[i]);
            supp.b[i] = std::min(supp_lambda.b[i],supp_mu.b[i]);
            if (supp.a[i] >= supp.b[i])
                return false;
        }
        return true;
    }


    template <class IBASIS, unsigned int DIM>
    void intersecting_wavelets(const TensorBasis<IBASIS,DIM>& basis,
                               const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                               const MultiIndex<int,DIM> j, const bool generators,
                               std::list<typename TensorBasis<IBASIS,DIM>::Index>& intersecting)
    {
        /*
         * we utilize neighboring relations of wavelets with nontrivial intersection in the translation index k!
         * 
         * iterate over all possible types on the level j. For all such j:
         * find the first valid k on this level and store its number -> current_num
         * find the first k with intersection -> k_hit
         * for all k's that differ only in k[dim] from k_hit:
         * - compute number of k's before k_hit ("jump_before[dim]" in the number)
         * - compute number of k's that intersect ("num_of_intersections[dim])
         * - compute number of k's after k_hit ("jump_after[dim]")
         * - add all wavelets from current_num+jump_before[dim] to current_num+jump_after[dim]
         * current_num = current_num + number of all k's that only differ in the last dimension
         * increase k_hit[dim-1] (if possible)
         * - if we have arrived at the last possible value for k_hit[dim-1]:
         * -- jump_after[dim-1] = number of possible values after k_hit[dim-1] times number of possible values for k[dim]
         * -- jump_before[dim-1] = number of possible values before k_hit[dim-1] times number of possible values for k[dim]
         * continue at num = currentnum + jump_after[dim-1]+jump_after[dim-1]
         * 
         * and analogous with dim-2 ...
         * this results in
         * 
         * while (!done)
         * {
         *      "jump_before" for all dimensions (current_jump_dim, current_jump_dim+1,...,dim)
         *      add wavelets
         *      increase k
         *      "jump_after" all dimensions where we are at the last possible value for k (dim, dim-1, ..., current_jump_dim)
         *      have we arrived at the very last possible k? -> increase type vector
         *      have we arrived at the type vector (1,1,...,1)? -> done = true
         * }
         * 
         * Also included: old code which works without get_intersecting_wavelets_on_level that is compatible with DSBasis
         * and even older code which works with brute force
         */
#if 1
        typedef typename TensorBasis<IBASIS,DIM>::Index Index;
        typedef typename IBASIS::Index Index1D;
        typedef typename Index::type_type type_type;
        if (generators) assert (j==basis.j0());
        intersecting.clear();
        // initialize the type-vector e
        type_type min_type, current_type;
        for (unsigned int i=0;i<DIM;i++)
        {
            min_type[i]=(j[i]==basis.j0()[i]) ?0:1; // remember whether we are on the minimal levels to avoid calls of basis.j0()
            current_type[i]=min_type[i];
        }
        if ( (!generators) && (j == basis.j0()))
        {
            current_type[DIM-1] = 1;
        }
        FixedArray1D<int,DIM> jump_before, jump_after, mink_gen, mink_wav, maxk_gen, maxk_wav;
        for (unsigned int i = 0; i < DIM; i++) 
        {
            if (min_type[i] == 0)
            {
                get_intersecting_wavelets_on_level(*(basis.bases()[i]),
                        lambda.j()[i],
                        lambda.e()[i],
                        lambda.k()[i],
                        j[i], true, mink_gen[i],maxk_gen[i]);
            }
            if (!(generators))
            {
                get_intersecting_wavelets_on_level(*basis.bases()[i],
                        lambda.j()[i],
                        lambda.e()[i],
                        lambda.k()[i],
                        j[i], false, mink_wav[i],maxk_wav[i]);
            }
        }
        unsigned int cjd(0); // current_jump_dim
        int blocksize;
        unsigned int number;
        if (generators == true)
        {
            number = 0;
        }
        else
        {
            number = basis.first_wavelet(j).number();
        }
        bool done = false;
        if (DIM == 1)
        {
            if (generators == true)
            {
                number += mink_gen[0] - basis.bases()[0]->DeltaLmin();
                blocksize = maxk_gen[0] - mink_gen[0] +1;
            } else
            {
                number += mink_wav[0] - basis.bases()[0]->Nablamin();
                blocksize = maxk_wav[0] - mink_wav[0] +1;
            }
    // add wavelets
            for (unsigned int n = number; n < number + blocksize; ++n)
            {
                intersecting.push_back(basis.get_wavelet(n));
            }
            done = true;
        }
        int basf(0); // "blocks added so far" (on the current type e) via modulo calculus we can deduce which block we have to add next
        while (!done)
        {
            if (cjd == 0)
            {
            // compute all jump_numbers that are needed for the currenttype
                int temp_int(1);
                for (int i = DIM-1; i >= 0; i--)
                {
                    if (current_type[i] == 0)
                    {
                        jump_before[i] = temp_int;
                        jump_before[i] *= mink_gen[i] - basis.bases()[i]->DeltaLmin();
                        jump_after[i] = temp_int;
                        jump_after[i] *= basis.bases()[i]->DeltaRmax(j[i]) - maxk_gen[i];
                        temp_int *= basis.bases()[i]->Deltasize(j[i]);
                    }
                    else
                    {
                        jump_before[i] = temp_int;
                        jump_before[i] *= mink_wav[i] - basis.bases()[i]->Nablamin();
                        jump_after[i] = temp_int;
                        jump_after[i] *= basis.bases()[i]->Nablamax(j[i]) - maxk_wav[i];
                        temp_int *= basis.bases()[i]->Nablasize(j[i]);
                    }
                }
                for (int i = DIM-2; i >= 0; --i)
                {
                    jump_before[i] += jump_before[i+1];
                    jump_after[i] += jump_after[i+1];
                }
                if (current_type[DIM-1] == 0)
                {
                    blocksize = maxk_gen[DIM-1] - mink_gen[DIM-1] +1;
                } else
                {
                    blocksize = maxk_wav[DIM-1] - mink_wav[DIM-1] +1;
                }
                basf = 0; // no blocks with the current type were added so far
            }
            // "jump_before" for all dimensions (current_jump_dim, current_jump_dim+1,...,dim)
            number += jump_before[cjd];
            // add wavelets
            for (unsigned int n = number; n < number + blocksize; ++n)
            {
                intersecting.push_back(basis.get_wavelet(n));
            }
            ++basf; // we have added a block!
            number += blocksize;
            // "increase k"
            // "jump_after" all dimensions where we are at the last possible value for k (dim, dim-1, ..., current_jump_dim)
            // have we arrived at the very last possible k? -> done = true
            
            // compute the dimension in which we want to increase k
            bool last_k(false); // have we just added the last possible block for the type vector e?
            int temp_i(basf);
            
            for (int i(DIM-1); i >= 0; --i)
            {
                if (i == 0) //we have already added the last possible block for this type, i.e., k[0] = maxk[0]
                {
                    last_k = true;
                    cjd = 0;
                    break;
                }
                int temp_count = (current_type[i-1] == 0)? (maxk_gen[i-1] - mink_gen[i-1] + 1) : (maxk_wav[i-1] - mink_wav[i-1] + 1); // blocksize in this dimension
                if ((temp_i % temp_count) != 0)
                {
                    cjd = i; // == "we can increase k[cjd]"
                    break;
                }
                temp_i = temp_i / temp_count;
            }
            // have we arrived at the very last possible k? -> increase type vector
            if (last_k == true)
            {
                done = true;
                if (generators == false) // == true => there is only one allowed type vector
                {
                    // "small loop"
                    // try to increase currenttype
                    // iterate over all combinations of generators/wavelets for all dimensions with currentlevel[i]=j0[i]
                    // this looks like binary addition: (in 3 Dim:) gwg is followed by gww (g=0,w=1)
                    for (int i(DIM-1); i >= 0; i--)
                    {
                        // find first position on level j0
                        if (j[i] == basis.j0()[i])
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
                }
            }
            // have we arrived at the type vector (1,1,...,1)? -> done = true
            // done == true means that all components with currentlevel[i]=j0[i] were wavelets.
            number += jump_after[cjd];
        }
#endif
#if 0 //compatible with DSBasis -> produces unsorted output
        
        if (generators) assert (j==basis.j0());

        typedef typename TensorBasis<IBASIS,DIM>::Index Index;
        intersecting.clear();

        // the set of intersecting wavelets is a cartesian product from d sets from the 1D case,
        // so we only have to compute the relevant 1D indices
        typedef typename IBASIS::Index Index1D;
        FixedArray1D<std::list<Index1D>,DIM> intersecting_1d_generators, intersecting_1d_wavelets;
        // prepare all intersecting wavelets and generators in the i-th coordinate direction
        for (unsigned int i = 0; i < DIM; i++) {
            if ((j[i]==basis.j0()[i]) && (generators || DIM > 1) )
                intersecting_wavelets(*basis.bases()[i],
                                      Index1D(lambda.j()[i],
                                              lambda.e()[i],
                                              lambda.k()[i],
                                              basis.bases()[i]),
                                      j[i], true, intersecting_1d_generators[i]);
            if (!generators)
                intersecting_wavelets(*basis.bases()[i],
                                      Index1D(lambda.j()[i],
                                              lambda.e()[i],
                                              lambda.k()[i],
                                              basis.bases()[i]),
                                      j[i], false, intersecting_1d_wavelets[i]);
        }
        // generate all relevant tensor product indices with either e=(0,...,0) or e!=(0,...,0)

        // Version 1: unsortierter Output
        // PERFORMANCE: Reihenfolge entspricht nicht der von tbasis_index. 
        // Problem? Ja! Reihenfolge wird in apply(window,y,res) ausgenutzt,
        // um heraus zu finden, ob der Zeileniterator erhöht werden muss
        typedef std::list<FixedArray1D<Index1D,DIM> > list_type;
        list_type indices;
        FixedArray1D<Index1D,DIM> helpindex;
        if ((j[0] == basis.j0()[0]) && (generators || DIM > 1) ) {
            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[0].begin()),
                    itend(intersecting_1d_generators[0].end());
                    it != itend; ++it) {
                        helpindex[0] = *it;
                        indices.push_back(helpindex);
                    }
        }
        if (!(generators)) {
            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[0].begin()),
                    itend(intersecting_1d_wavelets[0].end());
                    it != itend; ++it) {
                        helpindex[0] = *it;
                        indices.push_back(helpindex);
                    }
        }
        for (unsigned int i = 1; i < DIM; i++) {
            list_type sofar;
            sofar.swap(indices);
            for (typename list_type::const_iterator itready(sofar.begin()), itreadyend(sofar.end());
                    itready != itreadyend; ++itready) {
                        helpindex = *itready;
                        unsigned int esum = 0;
                        for (unsigned int k = 0; k < i; k++)
                            esum += helpindex[k].e();
                        if (generators || ( (j[i] == basis.j0()[i]) && ( i < DIM-1 || (i == (DIM-1) && esum > 0) ) ) ) {
                            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[i].begin()),
                                    itend(intersecting_1d_generators[i].end());
                                    it != itend; ++it) {
                                        helpindex[i] = *it;
                                        indices.push_back(helpindex);
                                    }
                        }
                        if (!(generators)) {
                            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[i].begin()),
                                    itend(intersecting_1d_wavelets[i].end());
                                    it != itend; ++it) {
                                        helpindex[i] = *it;
                                        indices.push_back(helpindex);
                                    }
                        }
                    }
        }
        // compose the results
        typename Index::type_type help_e;
        typename Index::translation_type help_k;
        for (typename list_type::const_iterator it(indices.begin()), itend(indices.end()); it != itend; ++it)
        {
            for (unsigned int i = 0; i < DIM; i++)
            {
                help_e[i] = (*it)[i].e();
                help_k[i] = (*it)[i].k();
            }
            intersecting.push_back(Index(j, help_e, help_k, &basis));
        }
#endif
        /*
        // a brute force solution
        typedef typename TensorBasis<IBASIS,DIM>::Support Support;
        Support supp;
        if (generators) {
            for (Index mu = first_generator<IBASIS,DIM>(&basis, j);; ++mu) {
                if (intersect_supports(basis, lambda, mu, supp))
                    intersecting.push_back(mu);
                if (mu == last_generator<IBASIS,DIM>(&basis, j)) break;
            }
        } else {
            for (Index mu = first_wavelet<IBASIS,DIM>(&basis, j);; ++mu) {
                if (intersect_supports(basis, lambda, mu, supp))
                    intersecting.push_back(mu);
                if (mu == last_wavelet<IBASIS,DIM>(&basis, j)) break;
            }
        }
         */
    }
    
#if 0
    template <class IBASIS, unsigned int DIM>
    void intersecting_elements(const TensorBasis<IBASIS,DIM>& basis,
                               const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                               const MultiIndex<int,DIM> j,
                               std::list<typename TensorBasis<IBASIS,DIM>::Index>& intersecting)
    {
        typedef typename TensorBasis<IBASIS,DIM>::Index Index;
        intersecting.clear();

        // the set of intersecting wavelets is a cartesian product from d sets from the 1D case,
        // so we only have to compute the relevant 1D indices
        typedef typename IBASIS::Index Index1D;
        FixedArray1D<std::list<Index1D>,DIM> intersecting_1d_generators, intersecting_1d_wavelets;
        // prepare all intersecting wavelets and generators in the i-th coordinate direction
        // initialize the type-vector e
        typename Index::type_type min_type,current_type;
        for (unsigned int i=0;i<DIM;i++)
        {
            min_type[i]=(j[i]==basis.j0()[i]) ?0:1;
            current_type[i]=min_type[i];
        }
        
        for (unsigned int i = 0; i < DIM; i++)
        {

            //if (mintype[i]==0)
            if (min_type[i] == 0)
                intersecting_wavelets(*basis.bases()[i],
                                      Index1D(lambda.j()[i],
                                              lambda.e()[i],
                                              lambda.k()[i],
                                              basis.bases()[i]),
                                      j[i], true, intersecting_1d_generators[i]);

            intersecting_wavelets(*basis.bases()[i],
                                  Index1D(lambda.j()[i],
                                          lambda.e()[i],
                                          lambda.k()[i],
                                          basis.bases()[i]),
                                  j[i], false, intersecting_1d_wavelets[i]);
        }

        // generate all tensor product indices
        // --------------------------------------------
        // speichere zu jedem Typ was bisher berechnet wurde. key = nummer des typs
        typedef std::list<FixedArray1D<Index1D,DIM> > list_type;
        typedef std::list<list_type> type_storage;
        type_storage storage;

        list_type temp_indices;
        FixedArray1D<Index1D,DIM> helpindex;
        if (min_type[0] == 0)
        {
            for (typename std::list<Index1D>::const_iterator it(intersecting_1d_generators[0].begin()),
                    itend(intersecting_1d_generators[0].end());
                    it != itend; ++it)
            {
                helpindex[0] = *it;
                temp_indices.push_back(helpindex);
            }
            storage.push_back(temp_indices);
            temp_indices.clear();
        }

        for (typename std::list<Index1D>::const_iterator it(intersecting_1d_wavelets[0].begin()),
                itend(intersecting_1d_wavelets[0].end());
                it != itend; ++it)
        {
            helpindex[0] = *it;
            temp_indices.push_back(helpindex);
        }
        storage.push_back(temp_indices);
        temp_indices.clear();

        for (int i = 1; i < DIM; i++)
        {
            type_storage sofar;
            sofar.swap(storage);
            for (typename type_storage::const_iterator it(sofar.begin()), itend(sofar.end()); it != itend; ++it)
            {
                // combine every element in "it" with every generator (if applicable) and wavelet in direction i
                //list_type sofar1;
                //sofar1.swap(indices1);
                if (min_type[i] == 0)
                { // add all possible generators
                    for (typename list_type::const_iterator itready((*it).begin()), itreadyend((*it).end()); itready != itreadyend; ++itready)
                    {
                        helpindex = *itready;
                        //unsigned int esum = 0;
                        //for (unsigned int k = 0; k < i; k++)
                        //esum += helpindex[k].e();
                        for (typename std::list<Index1D>::const_iterator it2(intersecting_1d_generators[i].begin()), itend2(intersecting_1d_generators[i].end()); it2 != itend2; ++it2)
                        {
                            helpindex[i] = *it2;
                            temp_indices.push_back(helpindex);
                        }
                    }
                    storage.push_back(temp_indices);
                    temp_indices.clear();
                }
                // now add all possible wavelets
                for (typename list_type::const_iterator itready((*it).begin()), itreadyend((*it).end()); itready != itreadyend; ++itready)
                {
                    helpindex = *itready;
                    for (typename std::list<Index1D>::const_iterator it2(intersecting_1d_wavelets[i].begin()), itend2(intersecting_1d_wavelets[i].end()); it2 != itend2; ++it2)
                    {
                        helpindex[i] = *it2;
                        temp_indices.push_back(helpindex);
                    }
                }
                storage.push_back(temp_indices);
                temp_indices.clear();
            }
        }
        // compose the results
        typename Index::type_type help_e;
        typename Index::translation_type help_k;
        for (typename type_storage::const_iterator it_store(storage.begin()),it_store_end(storage.end());it_store != it_store_end; it_store++)
        {
            //list_type somelist(*it_store);
            //typename list_type::const_iterator it22(*it_store.begin()), itend22(*it_store.end());
            for (typename list_type::const_iterator it((*it_store).begin()), itend((*it_store).end()); it != itend; ++it)
            {
                for (unsigned int i = 0; i < DIM; i++)
                {
                    help_e[i] = (*it)[i].e();
                    help_k[i] = (*it)[i].k();
                }
                intersecting.push_back(Index(j, help_e, help_k, &basis));
            }
        }
    }
#endif

    template <class IBASIS, unsigned int DIM>
    bool intersect_singular_support(const TensorBasis<IBASIS,DIM>& basis,
                                    const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                                    const typename TensorBasis<IBASIS,DIM>::Index& mu)
    {
        int j, k1, k2;
        for (unsigned int i = 0; i < DIM; i++) {
            if (!(intersect_singular_support(*basis.bases()[i],
                                             lambda.j()[i], lambda.e()[i], lambda.k()[i],
                                             mu.j()[i], mu.e()[i], mu.k()[i],
                                             j, k1, k2) ))
                return false;
        }
        return true;
    }
}

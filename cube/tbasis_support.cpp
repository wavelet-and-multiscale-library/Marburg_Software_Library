// implementation for tbasis_support.h

#include <utils/multiindex.h>
#include <utils/fixed_array1d.h>

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
        WaveletTL::support<IBASIS,DIM>(basis, lambda, supp_lambda);
        typename TensorBasis<IBASIS,DIM>::Support supp_mu;
        WaveletTL::support<IBASIS,DIM>(basis, mu, supp_mu);
        // determine support intersection granularity,
        // adjust single support granularities if necessary
        for (unsigned int i=0;i<DIM;i++)
        {
            supp.j[i] = std::max(supp_lambda.j[i], supp_mu.j[i]);
            if (supp_lambda.j[i] > supp_mu.j[i]) {
                const int adjust = 1<<(supp_lambda.j[i]-supp_mu.j[i]);
                supp_mu.a[i] *= adjust;
                supp_mu.b[i] *= adjust;
            } else {
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
        if (generators) assert (j==basis.j0());

        typedef typename TensorBasis<IBASIS,DIM>::Index Index;
        intersecting.clear();
#if 1
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
        // PERFORMANCE: Reihenfolge entspricht nicht der von tbasis_index. Problem?
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
        for (typename list_type::const_iterator it(indices.begin()), itend(indices.end());
                it != itend; ++it) {
                    for (unsigned int i = 0; i < DIM; i++) {
                        help_e[i] = (*it)[i].e();
                        help_k[i] = (*it)[i].k();
                    }
                    intersecting.push_back(Index(j, help_e, help_k, &basis));
                }
#else
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
#endif
    }

    template <class IBASIS, unsigned int DIM>
    bool intersect_singular_support(const TensorBasis<IBASIS,DIM>& basis,
                                    const typename TensorBasis<IBASIS,DIM>::Index& lambda,
                                    const typename TensorBasis<IBASIS,DIM>::Index& mu)
    {
        // we have intersection of the singular supports if and only if
        // (cube_support)   one of the components has this property in one dimension
        // (tbasis_support) all of the components have this property
        typedef typename IBASIS::Index Index1D;
        for (unsigned int i = 0; i < DIM; i++) {
            if (!(intersect_singular_support(*basis.bases()[i],
                                             Index1D(lambda.j()[i], lambda.e()[i], lambda.k()[i], basis.bases()[i]),
                                             Index1D(mu.j()[i], mu.e()[i], mu.k()[i], basis.bases()[i]))
                          ))
                return false;
        }
        return true;
    }
}


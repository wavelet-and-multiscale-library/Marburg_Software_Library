// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_ADAPTED_BASIS_H
#define _WAVELETTL_ADAPTED_BASIS_H

#include <interval/i_adapted_index.h>
#include <algebra/infinite_vector.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    Template wrapper class for adapting a multiwavelet basis as a classical basis
  */
  template <class IBASIS>
  class AdaptedBasis
  {
    public:
    /*!
      constructor with instance of basis to wrap (also serves as default constructor)
    */
    AdaptedBasis(IBASIS* basis = NULL)
      : basis_(basis)
    {
    }

    //! coarsest possible level
    inline const int j0() const { return basis_->j0(); }

    //! wavelet index class
    typedef IntervalAdaptedIndex<AdaptedBasis<IBASIS> > Index;

    //! unadapted index class
    typedef typename IBASIS::Index MultiIndex;

    //! size_type, for convenience
    typedef typename IBASIS::size_type size_type;

    //! geometric type of the support sets (2^{-j}[k1,k2])
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single primal generator or wavelet with index lambda.
      Note that wavelets have granularity j = lambda.j()+1, so the
      returned k values are respective to lambda.j()+lambda.e().
    */
    void primal_support(const Index& lambda, int& k1, int& k2) const
    { basis_->primal_support(*lambda.multi_index(), k1, k2); }

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single dual generator or wavelet with index lambda.
      Note that wavelets have granularity j = lambda.j()+1, so the
      returned k values are respective to lambda.j()+lambda.e().
    */
    void dual_support(const Index& lambda, int& k1, int& k2) const
    { basis_->dual_support(*lambda.multi_index(), k1, k2); }

    //! space dimension of the underlying domain
    static const int space_dimension = IBASIS::space_dimension;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return IBASIS::primal_regularity(); }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return IBASIS::primal_polynomial_degree(); }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }

    //! critical Sobolev regularity for the dual generators/wavelets
    static double dual_regularity() { return IBASIS::dual_regularity(); }

    //! degree of polynomial reproduction for the dual generators/wavelets
    static unsigned int dual_polynomial_degree() { return IBASIS::dual_polynomial_degree(); }

    //! number of vanishing moments for the dual wavelets
    static unsigned int dual_vanishing_moments() { return IBASIS::dual_vanishing_moments(); }


    /*!
      extremal generator indices in \Delta_j^X and \tilde\Delta_j^X
    */
    // no boundary indices on the primal side
    inline const int DeltaLmin() const { return IBASIS::Index::ck_encode(basis_->DeltaLmin(),0); }
    inline const int DeltaLmax() const { return IBASIS::Index::ck_encode(basis_->DeltaLmax(),IBASIS::number_of_components-1); }
    inline const int Delta0min() const { return DeltaLmax()+1; }
    inline const int Delta0max(const int j) const { return DeltaRmin(j)-1; }
    inline const int DeltaRmin(const int j) const { return IBASIS::Index::ck_encode(basis_->DeltaRmin(j),0); }
    inline const int DeltaRmax(const int j) const { return IBASIS::Index::ck_encode(basis_->DeltaRmax(j),IBASIS::number_of_components-1); }
    
    inline const int DeltaLTmin() const { return IBASIS::Index::ck_encode(basis_->DeltaLTmin(),0); } // == DeltaLmin()
    inline const int DeltaLTmax() const { return IBASIS::Index::ck_encode(basis_->DeltaLTmax(),IBASIS::number_of_components-1); }
    inline const int Delta0Tmin() const { return DeltaLTmax()+1; }
    inline const int Delta0Tmax(const int j) const { return DeltaRTmin(j)-1; }
    inline const int DeltaRTmin(const int j) const { return IBASIS::Index::ck_encode(basis_->DeltaRTmin(j),0); }
    inline const int DeltaRTmax(const int j) const { return IBASIS::Index::ck_encode(basis_->DeltaRmax(j),IBASIS::number_of_components-1); } // == DeltaRmax(j)

    //! size of Delta_j
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }
    
    /*!
      extremal wavelet indices in \nabla_j
    */
    inline const int Nablamin() const { return IBASIS::Index::ck_encode(basis_->Nablamin(),0); }
    inline const int Nablamax(const int j) const { return IBASIS::Index::ck_encode(basis_->Nablamax(j),IBASIS::number_of_components-1); }
    inline const int NablaLmax() const { return IBASIS::Index::ck_encode(basis_->NablaLmax(),IBASIS::number_of_components-1); }
    inline const int NablaRmin(const int j) const { return IBASIS::Index::ck_encode(basis_->NablaRmin(j),0); }

    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return Nablamax()-Nablamin()+1; }


    //! index of first (leftmost) generator on level j >= j0
    Index first_generator(const int j) const
    { return Index(basis_->first_generator(j), this); }

    //! index of last (rightmost) generator on level j >= j0
    Index last_generator(const int j) const
    { return Index(basis_->last_generator(j), this); }

    //! index of first (leftmost) wavelet on level j >= j0
    Index first_wavelet(const int j) const
    { return Index(basis_->first_wavelet(j), this); }

    //! index of last (rightmost) wavelet on level j >= j0
    Index last_wavelet(const int j) const
    { return Index(basis_->last_wavelet(j), this); }


    //! primal DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
      \psi_\lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where the multiscale decomposition starts with the coarsest
      generator level jmin, i.e. all generators have level jmin,
      wavelets may have higher levels.
    */
    void decompose_1(const Index& lambda, const int jmin,
                     InfiniteVector<double, Index>& c) const;

    //! primal DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one v
      such that
      \sum_{\lambda}c_\lambda\psi_\lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
      where generator level == jmin and wavelet level >= jmin
    */
    void decompose(const InfiniteVector<double, Index>& c, const int jmin,
                   InfiniteVector<double, Index>& v) const;

    //! dual DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
      \tilde\psi_\lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where the multiscale decomposition starts with the coarsest
      generator level jmin, i.e. all generators have level jmin,
      wavelets may have higher levels.
    */
    void decompose_t_1(const Index& lambda, const int jmin,
                       InfiniteVector<double, Index>& c) const;

    //! dual DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one v with wavelet level >= jmin,
      generator level == jmin, such that
      \sum_{\lambda}c_\lambda\tilde\psi_\lambda = \sum_{\lambda'}d_{\lambda'}\tilde\psi_{\lambda'}
    */
    void decompose_t(const InfiniteVector<double, Index>& c, const int jmin,
                     InfiniteVector<double, Index>& v) const;

    //! primal RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
      \psi_\lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct_1(const Index& lambda, const int j,
                       InfiniteVector<double, Index>& c) const;

    //! primal RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one v,
      such that
      \sum_{\lambda}c_\lambda\psi_\lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct(const InfiniteVector<double, Index>& c, const int j,
                     InfiniteVector<double, Index>& v) const;

    //! dual RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
      \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct_t_1(const Index& lambda, const int j,
                         InfiniteVector<double, Index>& c) const;

    //! dual RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one v,
      such that
      \sum_{\lambda}c_\lambda\tilde\psi_\lambda = \sum_{\lambda'}v_{\lambda'}\tilde\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct_t(const InfiniteVector<double, Index>& c, const int j,
                       InfiniteVector<double, Index>& v) const;


    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda
     */
    double primal_evaluate(const unsigned int derivative, const Index& lambda, const double x) const
    { return basis_->primal_evaluate(derivative, *lambda.multi_index(), x); }

    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void primal_evaluate(const unsigned int derivative, const Index& lambda,
                         const Array1D<double>& points, Array1D<double>& values) const
    { basis_->primal_evaluate(derivative, *lambda.multi_index(), points, values); }


    /*!
      reading access to underlying multiwavelet basis
     */
    const IBASIS* multi_basis() const { return basis_; }


  protected:
    //! hold an instance of the adapted multiwavelet basis
    const IBASIS* basis_;

    //! encode a coefficient set of multi-indices as a set of adapted indices
    void encode_infinite_vector(const InfiniteVector<double, MultiIndex>& c_multi, InfiniteVector<double, Index>& c) const;

    //! decode a coefficient set of adapted indices to a set of multi-indices
    void decode_infinite_vector(const InfiniteVector<double, Index>& c_adapt, InfiniteVector<double, MultiIndex>& c_multi) const;
  };
}

#include <interval/adapted_basis.cpp>

#endif // _WAVELETTL_ADAPTED_BASIS_H

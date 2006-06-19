// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_BASIS_H
#define _WAVELETTL_LDOMAIN_BASIS_H

#include <list>
#include <map>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <algebra/block_matrix.h>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>

#include <Ldomain/ldomain_index.h>

// for convenience, include also some functionality
#include <Ldomain/ldomain_support.h>
#include <Ldomain/ldomain_evaluate.h>

using std::list;
using MathTL::FixedArray1D;
using MathTL::InfiniteVector;
using MathTL::SparseMatrix;
using MathTL::BlockMatrix;

namespace WaveletTL
{
  /*!
    Template class for composite wavelet bases over the L-shaped domain
      (-1,1)^2 \ (0,1)^2
    such that the primal and the dual basis fulfill homogeneous Dirichlet b.c.'s
    on the outer domain boundary.

    References:
    [DS] Dahmen, Schneider:
         Composite Wavelet Bases for Operator Equations
	 Math. Comput. 68 (1999), 1533-1567
  */
  template <class IBASIS>
  class LDomainBasis
  {
  public:
    //! default constructor
    LDomainBasis();

    //! interval basis
    typedef IBASIS IntervalBasis;

    //! coarsest possible level j0
    inline const int j0() const { return basis1d_.j0(); }
    
    //! wavelet index class
    typedef LDomainIndex<IBASIS> Index;

    //! geometric type of the support sets
    typedef struct {
      int j;       // granularity
      int xmin[3];
      int xmax[3];
      int ymin[3];
      int ymax[3];
    } Support;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return IBASIS::primal_regularity(); } // dirty, we should use max(1.5,~) instead
    
    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }

    //! number of vanishing moments for the dual wavelets
    static unsigned int dual_vanishing_moments() { return IBASIS::dual_vanishing_moments(); }

    //! read access to the underlying 1D basis
    const IntervalBasis& basis1d() const { return basis1d_; }

    //! size of Delta_j
    const int Deltasize(const int j) const;

    //! sizes of the different wavelet index sets
    const int Nabla01size(const int j) const;
    const int Nabla10size(const int j) const;
    const int Nabla11size(const int j) const;

    /*!
      The following routines provide read access to the diverse refinement (sub)matrices
      on a level j >= j0. The row and column indices follow the less<Index> ordering.
      Those matrices will be collected in an internal cache to provide faster access.
    */
    const BlockMatrix<double>&  get_Mj0      (const int j) const;
    const BlockMatrix<double>&  get_Mj0T     (const int j) const;
    const SparseMatrix<double>& get_Mj1c_1d  (const int j) const; // initial stable completion in 1D
    const BlockMatrix<double>&  get_Mj1c_01  (const int j) const;
    const BlockMatrix<double>&  get_Mj1c_10  (const int j) const;
    const BlockMatrix<double>&  get_Mj1c_11  (const int j) const;

    //! index of first generator on level j >= j0
    Index first_generator(const int j) const;
      
    //! index of last generator on level j >= j0
    Index last_generator(const int j) const;
      
    //! index of first wavelet on level j >= j0
    Index first_wavelet(const int j) const;
      
    //! index of first wavelet on level j >= j0 with type e
    Index first_wavelet(const int j, const typename Index::type_type& e) const;

    //! index of last wavelet on level j >= j0
    Index last_wavelet(const int j) const;

    //! RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
      \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct_1(const Index& lambda, const int j,
		       InfiniteVector<double, Index>& c) const;

    //! RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one v,
      such that
      \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
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

  protected:
    //! the interval 1d wavelet basis
    IntervalBasis basis1d_;

    //! caches for the diverse refinement matrices
    typedef std::map<int,BlockMatrix<double> > MatrixCache;
    mutable MatrixCache Mj0_cache, Mj0T_cache,
      Mj1c_01_cache, Mj1c_10_cache, Mj1c_11_cache;

    typedef std::map<int,SparseMatrix<double> > Matrix1DCache;
    mutable Matrix1DCache Mj1c_1d_cache;
  };
}

#include <Ldomain/ldomain_basis.cpp>

#endif

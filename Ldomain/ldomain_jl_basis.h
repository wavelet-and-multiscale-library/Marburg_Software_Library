// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_JL_BASIS_H
#define _WAVELETTL_LDOMAIN_JL_BASIS_H

#include <Ldomain/ldomain_jl_index.h>
#include <Rd/r_mw_index.h>
#include <algebra/infinite_vector.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    A wavelet basis on the L-shaped domain
      (-1,1)^2 \ (0,1)^2
    with homogeneous boundary conditions. The generators are given as
    cubic Hermite interpolatory splines, and the wavelets are essentially taken
    from the multiwavelet basis of [JL].
  */
  class LDomainJLBasis
  {
  public:
    //! default constructor
    LDomainJLBasis();

    //! coarsest possible level j0
    inline const int j0() const { return j0_; }
    
    //! wavelet index class
    typedef LDomainJLIndex Index;

    /*!
      geometric type of the support sets;
      all supports and intersections thereof can be written as the intersection
      of a rectangle with the L-shaped domain.
    */
    typedef struct {
      int j;       // granularity
      int xmin;
      int xmax;
      int ymin;
      int ymax;
    } Support;

    //! space dimension of the underlying domain
    static const int space_dimension = 2;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return 2.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return 4; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return 2; }

    //! RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
      \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct_1(const Index& lambda, const int j,
		       InfiniteVector<double,Index>& c) const;

    /*!
      helper function for 1D reconstruct_1() calls
    */
    void reconstruct_1_1d(const RMWIndex& lambda, const int j,
			  InfiniteVector<double,RMWIndex>& c) const;

//     //! RECONSTRUCT routine, full version
//     /*!
//       Constructs for a given coefficient set c another one v,
//       such that
//       \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
//       where always |\lambda'|>=j
//     */
//     void reconstruct(const InfiniteVector<double, Index>& c, const int j,
// 		     InfiniteVector<double, Index>& v) const;

    //! index of first generator on level j >= j0
    Index first_generator(const int j) const;
      
    //! index of last generator on level j >= j0
    Index last_generator(const int j) const;
      
    //! index of first wavelet on level j >= j0
    Index first_wavelet(const int j) const;
      
    //! index of first wavelet on level j >= j0 with type e
    Index first_wavelet(const int j, const Index::type_type& e) const;

    //! index of last wavelet on level j >= j0
    Index last_wavelet(const int j) const;

  protected:
    //! coarsest level
    int j0_;
  };
}

#include <Ldomain/ldomain_jl_basis.cpp>

#endif

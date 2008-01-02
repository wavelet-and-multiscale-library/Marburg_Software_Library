// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_JL_BASIS_H
#define _WAVELETTL_JL_BASIS_H

#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <interval/jl_index.h>
#include <utils/array1d.h>

using MathTL::Vector;
using MathTL::Matrix;
using MathTL::InfiniteVector;
using MathTL::Array1D;

namespace WaveletTL
{
  /*!
    Template class for the (multi)wavelet bases on the interval as introduced in [JL].
    Essentially, the primal generators

      phi_{j,0,1},...,phi_{j,0,2^j-1}  <-> 2^{j/2}phi_0(2^j*x-k), k=1,...,2^j-1
      phi_{j,1,0},...,phi_{j,0,2^j}    <-> 2^{j/2}phi_1(2^j*x-k), k=0,...,2^j

    are dilated and translated versions of the two cubic Hermite interpolants
    at 2^{-j}k, i.e.,
    
      phi_0(k) = delta_{0,k}, (d/dx)phi_0(k) = 0
      phi_1(k) = 0          , (d/dx)phi_1(k) = delta_{0,k}

    The wavelets are

      psi_{j,0,1},...,psi_{j,0,2^j-1}  <-> 2^{j/2}psi_0(2^j*x-k), k=1,...,2^j-1
      psi_{j,1,0},...,psi_{j,1,2^j}    <-> 2^{j/2}psi_1(2^j*x-k), k=0,...,2^j

    see [JL] for their definition.

    Note that here, in contrast to [JL], both components of the generators and
    the wavelets are normalized in the L_2 norm.

    References:
    [JL] Jia/Liu:
         Wavelet bases of Hermite cubic splines on the interval
  */
  class JLBasis
  {
  public:
    /*!
      constructor
    */
    JLBasis();

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef JLIndex Index;

    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! geometric type of the support sets
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    //! space dimension of the underlying domain
    static const int space_dimension = 1;

    void set_jmax(const int jmax) {
      jmax_ = jmax;
      setup_full_collection();
    }

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return 2.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return 4; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return 2; }

    //! size of Delta_j
    inline const int Deltasize(const int j) const { return 1<<(j+1); }

    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return 1<<(j+1); }
    
    //! index of first (leftmost) generator on level j >= j0
    Index first_generator(const int j) const;

    //! index of last (rightmost) generator on level j >= j0
    Index last_generator(const int j) const;

    //! index of first (leftmost) wavelet on level j >= j0
    Index first_wavelet(const int j) const;

    //! index of last (rightmost) wavelet on level j >= j0
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

    //! get the wavelet index corresponding to a specified number
    const inline Index* get_wavelet (const int number) const {
      return &full_collection[number];
    }

    //! number of wavelets between coarsest and finest level
    const int degrees_of_freedom() const { return full_collection.size(); };


  protected:
    //! coarsest possible level
    int j0_;

    //! finest possible level
    int jmax_;

    //! boundary condition orders at 0 and 1
    int s0, s1;

    //! general setup routine which is shared by the different constructors
    void setup();

    //! setup full collectin of wavelets between j0_ and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();

    //! collection of all wavelets between coarsest and finest level
    Array1D<Index> full_collection;


  };
}

#include <interval/jl_basis.cpp>

#endif

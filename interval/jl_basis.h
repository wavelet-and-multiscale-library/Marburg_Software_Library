// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_JL_BASIS_H
#define _WAVELETTL_JL_BASIS_H

#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>
#include <interval/i_index.h>

// for convenience, include also some functionality
// #include <interval/jl_evaluate.h>

using MathTL::Vector;
using MathTL::Matrix;

namespace WaveletTL
{
  /*!
    Template class for the (multi)wavelet bases on the interval as introduced in [JL].
    Essentially, the primal generators

      phi_{j,s0},...,phi_{j,2^j-s1}             <-> 2^{j/2}phi_0(2^j*x-k), k=s0,...,s^j-s1
      phi_{j,2^j-s1+1},...,phi_{j,2^{j+1}-s1+1} <-> 2^{j/2}phi_1(2^j*x-k), k=0,...,2^j

    are dilated and translated versions of the two cubic Hermite interpolants
    at 2^{-j}k, i.e.,

      phi_0(k) = delta_{0,k}, (d/dx)phi_0(k) = 0
      phi_1(k) = 0          , (d/dx)phi_1(k) = delta_{0,k}

    The wavelets are

      psi_{j,1},...,psi_{j,2^j-1}     <-> 2^{j/2}psi_0(2^j*x-k), k=1,...,2^j-1
      psi_{j,2^j},...,psi_{j,2^{j+1}} <-> 2^{j/2}psi_1(2^j*x-k), k=0,...,2^j

    see [JL] for their definition.

    References:
    [JL] Jia/Liu:
         Wavelet bases of Hermite cubic splines on the interval
  */
  class JLBasis
  {
  public:
    /*!
      constructor
      
      At the moment, you may (only) specify the order of the primal (s) boundary conditions
      at the left and right end of the interval [0,1].
    */
    JLBasis(const int s0 = 1, const int s1 = 1);

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef IntervalIndex<JLBasis> Index;

    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! geometric type of the support sets
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    //! read access to the primal b.c. order at x=0
    const int get_s0() const { return s0; }

    //! read access to the primal b.c. order at x=1
    const int get_s1() const { return s1; }

    //! extremal generator indices
    inline const int DeltaLmin() const { return s0; }
    inline const int DeltaRmax(const int j) const { return (1<<(j+1))-s1+1; }
    
    //! size of Delta_j
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }

    //! boundary indices in \nabla_j
    inline const int Nablamin() const { return 1; }
    inline const int Nablamax(const int j) const { return 1<<(j+1); }

    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return 1<<(j+1); }

  protected:
    //! coarsest possible level
    int j0_;

    //! boundary condition orders at 0 and 1
    int s0, s1;

    //! general setup routine which is shared by the different constructors
    void setup();
  };
}

#include <interval/jl_basis.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_S_BASIS_H
#define _WAVELETTL_S_BASIS_H

#include <utils/array1d.h>
#include <algebra/vector.h>
#include <algebra/matrix.h>
#include <algebra/infinite_vector.h>
#include <algebra/fixed_matrix.h>
#include <interval/i_multi_index.h>
#include <Rd/dhjk_mask.h>

using MathTL::Array1D;
using MathTL::Vector;
using MathTL::Matrix;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    Template class for the Hermite spline (multi)wavelet bases on the interval.

    Essentially, the primal generators

      phi_{j,i,1},...,phi_{j,i,2^j-1} <-> 2^{j/2} phi_i(2^j*x-k), k=1,...,2^j-1

    are dilated and translated versions of the two cubic Hermite interpolants
    at 2^{-j}k, i.e.,
    
      phi_0(k) = delta_{0,k}, (d/dx)phi_0(k) = 0
      phi_1(k) = 0          , (d/dx)phi_1(k) = delta_{0,k}

    References:
    [S] A. Schneider: Konstruktion von Multiwavelets und Anwendungen bei
        adaptiven numerischen Verfahren, Diplomarbeit, 2007
  */
  class SBasis
  {
  public:
    /*!
      constructor
    */
    SBasis();

    /*!
      constructor with boundary (primal) condition order as parameter
      present just for the signature, needed by CubeBasis
      Requires s0 and s1 to be 2 as this basis has fixed boundary conditions.
    */
    SBasis(const int s0, const int s1);

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef IntervalMultiIndex< SBasis > Index;

    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

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
    void primal_support(const Index& lambda, int& k1, int& k2) const;

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single dual generator or wavelet with index lambda.
      Note that wavelets have granularity j = lambda.j()+1, so the
      returned k values are respective to lambda.j()+lambda.e().
    */
    void dual_support(const Index& lambda, int& k1, int& k2) const;

    //! space dimension of the underlying domain
    static const int space_dimension = 1;

    //! number of compononents of the multi-wavelet basis
    static const unsigned int number_of_components = 2;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return 2.5; }

    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return 4; }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return 2; }

    //! critical Sobolev regularity for the dual generators/wavelets
    static double dual_regularity() { return 0.824926; }

    //! degree of polynomial reproduction for the dual generators/wavelets
    static unsigned int dual_polynomial_degree() { return 2; }

    //! number of vanishing moments for the dual wavelets
    static unsigned int dual_vanishing_moments() { return 4; }

    //! read access to the primal b.c. order at x=0
    const int get_s0() const { return 2; }

    //! read access to the primal b.c. order at x=1
    const int get_s1() const { return 2; }

    //! read access to the dual b.c. order at x=0
    const int get_sT0() const { return 0; }

    //! read access to the dual b.c. order at x=1
    const int get_sT1() const { return 0; }


    /*!
      extremal generator indices in \Delta_j^X and \tilde\Delta_j^X
    */
    // no boundary indices on the primal side
    inline const int DeltaLmin() const { return 1; }
    inline const int DeltaLmax() const { return 1; }
    inline const int Delta0min() const { return 1; }
    inline const int Delta0max(const int j) const { return (1<<j)-1; }
    inline const int DeltaRmin(const int j) const { return (1<<j)-1; }
    inline const int DeltaRmax(const int j) const { return (1<<j)-1; }
    
    inline const int DeltaLTmin() const { return 1; } // == DeltaLmin()
    inline const int DeltaLTmax() const { return 1; }
    inline const int Delta0Tmin() const { return DeltaLTmax()+1; }
    inline const int Delta0Tmax(const int j) const { return DeltaRTmin(j)-1; }
    inline const int DeltaRTmin(const int j) const { return (1<<j)-1; }
    inline const int DeltaRTmax(const int j) const { return (1<<j)-1; } // == DeltaRmax()

    //! size of Delta_j
    inline const int Deltasize(const int j) const { return DeltaRmax(j)-DeltaLmin()+1; }
    
    /*!
      extremal wavelet indices in \nabla_j
    */
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }
    inline const int NablaLmax() const { return Nablamin(); }
    inline const int NablaRmin(const int j) const { return Nablamax(j); }

    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return 1<<j; }


    //! index of first (leftmost) generator on level j >= j0
    Index first_generator(const int j) const;

    //! index of last (rightmost) generator on level j >= j0
    Index last_generator(const int j) const;

    //! index of first (leftmost) wavelet on level j >= j0
    Index first_wavelet(const int j) const;

    //! index of last (rightmost) wavelet on level j >= j0
    Index last_wavelet(const int j) const;


    //! primal DECOMPOSE routine, simple version, one step
    /*!
      Constructs for a given single generator index lambda a coefficient set c,
      such that
      \phi_\lambda = \sum_{\lambda'}c_{\lambda'}\phi_{\lambda'} + \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where |\lambda'| = |\lambda|-1
    */
    void decompose_1(const Index& lambda, InfiniteVector<double, Index>& c) const;

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

    //! dual DECOMPOSE routine, simple version, one step
    /*!
      Constructs for a given single generator index lambda a coefficient set c,
      such that
      \tilde\phi_\lambda = \sum_{\lambda'}c_{\lambda'}\tilde\phi_{\lambda'} + \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where |\lambda'| = |\lambda|-1
    */
    void decompose_t_1(const Index& lambda, InfiniteVector<double, Index>& c) const;

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

    //! primal RECONSTRUCT routine, simple version, one step
    /*!
       for a given single wavelet index, write it as a sum of generators
        \psi_\lambda = \sum_{\lambda'} c_{\lambda'} \phi_{\lambda'}
       with |\lambda'| = |\lambda|+1
     */
    void reconstruct_1(const Index& lambda, InfiniteVector<double, Index>& c) const;

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

    //! dual RECONSTRUCT routine, simple version, one step
    /*!
       for a given single wavelet index, write it as a sum of generators
        \tilde\psi_lambda = \sum_{\lambda'} c_{\lambda'} \tilde\phi_{\lambda'}
       with |\lambda'| = |\lambda|+1
     */
    void reconstruct_t_1(const Index& lambda, InfiniteVector<double, Index>& c) const;

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

    //! primal REFINE routine, simple version, one step
    /*!
      For a given single primal generator index, compute a coefficient set c,
      such that
       \phi_lambda = \sum_{\lambda'} c_{\lambda'} \phi_{\lambda'}
      where |\lambda'| = |\lambda|+1
    */
    void refine_1(const Index& lambda, InfiniteVector<double, Index>& c) const;

    //! primal REFINE routine, simple version, until granularity j
    /*!
      For a given single primal generator index, compute a coefficient set c,
      such that
       \phi_lambda = \sum_{\lambda'} c_{\lambda'} \phi_{\lambda'}
      where always |\lambda'| >= j
    */
    void refine_1(const Index& lambda, const int j, InfiniteVector<double, Index>& c) const;

    //! primal REFINE routine, full version, until granularity j
    /*!
      For a given coefficient set c, compute another one v,
      such that
       \sum_{\lambda} c_\lambda \phi_lambda = \sum_{\lambda'} v_{\lambda'} \phi_{\lambda'}
      where always |\lambda'| >= j
    */
    void refine(const InfiniteVector<double, Index>& c, const int j,
                InfiniteVector<double, Index>& v) const;

    //! dual REFINE routine, simple version, one step
    /*!
      For a given single dual generator index, compute a coefficient set c,
      such that
       \tilde\phi_lambda = \sum_{\lambda'} c_{\lambda'} \tilde\phi_{\lambda'}
      where |\lambda'| = |\lambda|+1
    */
    void refine_t_1(const Index& lambda, InfiniteVector<double, Index>& c) const;

    //! dual REFINE routine, simple version, until granularity j
    /*!
      For a given single dual generator index, compute a coefficient set c,
      such that
       \tilde\phi_lambda = \sum_{\lambda'} c_{\lambda'} \tilde\phi_{\lambda'}
      where always |\lambda'| >= j
    */
    void refine_t_1(const Index& lambda, const int j, InfiniteVector<double, Index>& c) const;

    //! dual REFINE routine, full version, until granularity j
    /*!
      For a given coefficient set c, compute another one v,
      such that
       \sum_{\lambda} c_\lambda \tilde\phi_lambda = \sum_{\lambda'} v_{\lambda'} \tilde\phi_{\lambda'}
      where always |\lambda'| >= j
    */
    void refine_t(const InfiniteVector<double, Index>& c, const int j,
                  InfiniteVector<double, Index>& v) const;


    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda
     */
    double primal_evaluate(const unsigned int derivative, const Index& lambda, const double x) const;

    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    void primal_evaluate(const unsigned int derivative, const Index& lambda,
                         const Array1D<double>& points, Array1D<double>& values) const;


  protected:
    //! coarsest possible level
    static const int j0_ = 3;

    //! refinement masks
    DHJKMask_primal mask_primal;
    DHJKMask_dual mask_dual;

    //! dual refinement matrix boundary blocks MLT, MRT
    FixedMatrix<double, 8, 2> MLT; // left boundary
    FixedMatrix<double, 8, 2> MRT; // right boundary

    //! stable completion matrix blocks
    FixedMatrix<double, 6, 2> Mj1L; // left boundary
    FixedMatrix<double, 6, 2> Mj1R; // right boundary
    FixedMatrix<double, 10, 2> Mj1I; // inner block
    FixedMatrix<double, 6, 2> MTj1I; // has only inner blocks

    //! matrix block of the inverse
    FixedMatrix<double, 4, 2> GW;

    //! general setup routine which is shared by the different constructors
    void setup();
  };
}

#include <interval/s_basis.cpp>

#endif  // _WAVELETTL_S_BASIS_H

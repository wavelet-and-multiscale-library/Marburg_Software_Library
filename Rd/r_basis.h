// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_R_BASIS_H
#define _WAVELETTL_R_BASIS_H

#include <Rd/refinable.h>
#include <Rd/r_index.h>
#include <algebra/infinite_vector.h>
#include <utils/tiny_tools.h>

namespace WaveletTL
{
  /*!
    Base class for (bi)orthogonal wavelet bases in L_2(\mathbb R).

    Minimal signature of the mask template parameters modelling 1D refinable functions:
    
    const int abegin() const;          // inherited from RRefinementMask
    const int aend() const;            // dito
    const double a(const int k) const; // dito
    static double regularity() const;
    static unsigned int primal_polynomial_order();
  */
  template <class PRIMALMASK, class DUALMASK = PRIMALMASK>
  class RBasis
  {
  public:
    //! default constructor
    RBasis();

    //! wavelet index class
    typedef RIndex Index;

    //! make template arguments accessible
    typedef PRIMALMASK primal_mask;
    typedef DUALMASK   dual_mask;

    /*!
      coarsest possible level j0: we restrict ourselves to j>=0
    */
    static const int j0() { return 0; }

    //! critical Sobolev regularity for the primal generators/wavelets
    static inline double primal_regularity() { return PRIMALMASK::regularity(); }

    /*/
      degree of polynomial reproduction for the primal generators/wavelets
      == number of vanishing moments for the dual wavelets
    */
    static unsigned int primal_polynomial_degree() { return DUALMASK::Strang_Fix_order(); }

    //! number of vanishing moments for the primal wavelets
    static inline unsigned int primal_vanishing_moments() { return DUALMASK::Strang_Fix_order(); }

    //! reading access to the primal refinement mask (a_k)
    inline const PRIMALMASK& a() const { return primal_mask_; }

    //! reading access to the primal refinement mask (a_k)
    inline const double a(const int k) const { return a().a(k); }

    //! get primal mask (a_k) and start offset
    void get_a(Array1D<double>& a, int& offset) const;

    //! reading access to the dual refinement mask (aT_k)
    inline const DUALMASK& aT() const { return dual_mask_; }

    //! reading access to the dual refinement mask (aT_k)
    inline const double aT(const int k) const { return aT().a(k); }

    //! get dual mask (aT_k) and start offset
    void get_aT(Array1D<double>& aT, int& offset) const;

    //! reading access to the primal wavelet mask (b_k)
    inline const double b(const int k) const { return minus1power(k)*aT(1-k); }

    //! get primal wavelet mask (b_k) and start offset
    void get_b(Array1D<double>& b, int& offset) const;

    //! reading access to the dual wavelet coefficients
    inline const double bT(const int k) const { return minus1power(k)*a(1-k); }
    
    //! get dual wavelet mask (bT_k) and start offset
    void get_bT(Array1D<double>& bT, int& offset) const;
    
    /*!
      geometric type of the support sets
      (note that if k1 > k2, the support consists of two intervals instead of one)
    */
    typedef struct {
      int j;
      int k1;
      int k2;
    } Support;

    /*!
      Compute an interval 2^{-j}[k1,k2] which contains the support of a
      single primal spline generator or wavelet \psi_\lambda.
      (j == lambda.j()+lambda.e() is implicitly assumed for performance reasons)
    */
    static void support(const Index& lambda, int& k1, int& k2);

    //! DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where the multiscale decomposition starts with coarsest
      generator level j0. Each concrete instance of RBasis has to
      implement this routine.
     */
    void decompose_1(const Index& lambda, const int j0,
		     InfiniteVector<double, Index>& c) const;

    //! dual DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where the multiscale decomposition starts with coarsest
      generator level j0. Each concrete instance of RBasis has to
      implement this routine.
      (note that we don't use the method name DECOMPOSEt, which would irritate
      the linker)
     */
    void decompose_t_1(const Index& lambda, const int j0,
		       InfiniteVector<double, Index>& c) const;

    //! DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one d with level >= j0,
      such that
        \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}d_{\lambda'}\psi_{\lambda'}
    */
    void decompose(const InfiniteVector<double, Index>& c, const int j0,
		   InfiniteVector<double, Index>& d) const;

    //! dual DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one d with level >= j0,
      such that
        \sum_{\lambda}c_\lambda\tilde\psi_lambda = \sum_{\lambda'}d_{\lambda'}\tilde\psi_{\lambda'}
    */
    void decompose_t(const InfiniteVector<double, Index>& c, const int j0,
		     InfiniteVector<double, Index>& d) const;

    //! RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
      Each concrete instance of RBasis has to implement this routine.
      (note that we again don't use the method name RECONSTRUCT, which would irritate
      the linker)
     */
    void reconstruct_1(const Index& lambda, const int j,
		       InfiniteVector<double, Index>& c) const;

    //! RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one d,
      such that
        \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}d_{\lambda'}\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct(const InfiniteVector<double, Index>& c, const int j,
		     InfiniteVector<double, Index>& d) const;

    //! dual RECONSTRUCT routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where always |\lambda'|>=j
      Each concrete instance of RBasis has to implement this routine.
      (note that we again don't use the method name RECONSTRUCT, which would irritate
      the linker)
     */
    void reconstruct_t_1(const Index& lambda, const int j,
			 InfiniteVector<double, Index>& c) const;

    //! dual RECONSTRUCT routine, full version
    /*!
      Constructs for a given coefficient set c another one d,
      such that
        \sum_{\lambda}c_\lambda\tilde\psi_\lambda = \sum_{\lambda'}d_{\lambda'}\tilde\psi_{\lambda'}
      where always |\lambda'|>=j
    */
    void reconstruct_t(const InfiniteVector<double, Index>& c, const int j,
		       InfiniteVector<double, Index>& d) const;

    //! evaluate n-th derivative of a single wavelet on a dyadic grid
    /*!
      Evaluate the n-th derivative of a single primal/dual generator or wavelet \psi_\lambda
      on a dyadic subgrid of the interval [A,B].
      We assume that the derivative vanishes at the boundary of its support.
     */
    SampledMapping<1> evaluate(const unsigned int derivative,
			       const Index& lambda,
			       const bool primal,
			       const int A, const int B,
			       const int resolution) const;
    
    //! evaluate n-th derivative of a linear combination of wavelets
    /*!
      Evaluate an arbitrary linear combination of primal or dual
      wavelets on a dyadic subgrid of [A,B].
    */
    SampledMapping<1> evaluate(const unsigned int derivative,
			       const InfiniteVector<double, Index>& coeffs,
			       const bool primal,
			       const int A, const int B,
			       const int resolution) const;

    //! bands and offsets for the 4 masks
    int offset_a, offset_aT, offset_b, offset_bT;
    Array1D<double> band_a, band_aT, band_b, band_bT;
    
  protected:
    //! one instance of the primal mask
    PRIMALMASK primal_mask_;
    
    //! one instance of the dual mask
    DUALMASK dual_mask_;
  };
}

// include implementation
#include <Rd/r_basis.cpp>

#endif

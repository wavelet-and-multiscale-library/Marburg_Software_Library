// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_R_BASIS_H
#define _WAVELETTL_R_BASIS_H

#include <Rd/refinable.h>
#include <Rd/r_index.h>
#include <algebra/infinite_vector.h>

namespace WaveletTL
{
  /*!
    Abstract base class for (bi)orthogonal wavelet bases in L_2(\mathbb R).
  */
  template <class PRIMALMASK, class DUALMASK = PRIMALMASK>
  class RBasis
  {
  public:
    /*!
      default constructor
    */
    RBasis();

    /*!
      wavelet index class
    */
    typedef RIndex Index;

    /*!
      reading access to the primal mask
    */
    inline const RefinableFunction<PRIMALMASK> a() const { return a_; }

    /*!
      reading access to the dual mask
    */
    inline const RefinableFunction<DUALMASK> aT() const { return aT_; }

    /*!
      reading access to the primal wavelet coefficients
    */
    inline const MultivariateLaurentPolynomial<double, 1>& b() const { return b_; }

    /*!
      reading access to the dual wavelet coefficients
    */
    inline const MultivariateLaurentPolynomial<double, 1>& bT() const { return bT_; }

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
      Constructs for a given single wavelet index lambda a coefficient set d,
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
      Constructs for a given single wavelet index lambda a coefficient set d,
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

  protected:
    /*!
      one instance of the primal mask
    */
    RefinableFunction<PRIMALMASK> a_;

    /*!
      one instance of the dual mask
    */
    RefinableFunction<DUALMASK> aT_;

    /*!
      primal wavelet coefficients
    */
    MultivariateLaurentPolynomial<double, 1> b_;

    /*!
      dual wavelet coefficients
    */
    MultivariateLaurentPolynomial<double, 1> bT_;
  };
}

// include implementation
#include <Rd/r_basis.cpp>

#endif

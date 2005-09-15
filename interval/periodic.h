// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PERIODIC_H
#define _WAVELETTL_PERIODIC_H

#include <iostream>

#include <algebra/infinite_vector.h>

#include <interval/i_index.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    Template class for a periodic, biorthogonal wavelet basis on the unit interval [0,1].
    The template parameters provide the biorthogonal masks of a refinable function
    on the real line.
    If the second template parameter is omitted, we assume that the given refinable
    function induces an orthonormal wavelet basis (i.e. a biorthogonal one where the dual
    basis is just the primal one).

    The periodized scaling functions (and analogously the wavelets) look like

      phi^per_{j,k}(x) = \sum_{l\in\mathbb Z} phi_{j,k}(x+l)

    A univariate mask should have at least the signature of a LaurentPolynomial<double>,
    i.e., have the usual iterator classes. Typical examples are the Haar or the CDF masks.

    References:
    [D] Daubechies,
        Ten Lectures On Wavelets, pp. 304ff
  */
  template <class PRIMALMASK, class DUALMASK = PRIMALMASK>
  class PeriodicBasis
  {
  public:
    //! default constructor
    PeriodicBasis();
    
    //! coarsest possible level, this is always zero
    inline const int j0() const { return 0; }

    //! wavelet index class
    typedef IIndex<PeriodicBasis<PRIMALMASK,DUALMASK> > Index;

    //! bounds for the generator indices
    inline const int DeltaLmin() const { return 0; }
    inline const int DeltaRmax(const int j) const { return (1<<j) - 1; }

    //! bounds for the wavelet indices
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j) - 1; }
    
    //! size of Delta_j
    inline const int Deltasize(const int j) const { 1<<j; }
    
    //! DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \psi_lambda = \sum_{\lambda'}c_{\lambda'}\psi_{\lambda'}
      where the multiscale decomposition starts with the coarsest
      generator level jmin.
     */
    void decompose_1(const Index& lambda, const int jmin,
		     InfiniteVector<double, Index>& c) const;

    //! dual DECOMPOSE routine, simple version
    /*!
      Constructs for a given single wavelet index lambda a coefficient set c,
      such that
        \tilde\psi_lambda = \sum_{\lambda'}c_{\lambda'}\tilde\psi_{\lambda'}
      where the multiscale decomposition starts with the coarsest
      generator level jmin.
     */
    void decompose_t_1(const Index& lambda, const int jmin,
		       InfiniteVector<double, Index>& c) const;

    //! DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one v with level >= jmin,
      such that
        \sum_{\lambda}c_\lambda\psi_lambda = \sum_{\lambda'}v_{\lambda'}\psi_{\lambda'}
    */
    void decompose(const InfiniteVector<double, Index>& c, const int jmin,
		   InfiniteVector<double, Index>& v) const;

    //! dual DECOMPOSE routine, full version
    /*!
      constructs for a given coefficient set c another one v with level >= jmin,
      such that
        \sum_{\lambda}c_\lambda\tilde\psi_lambda = \sum_{\lambda'}d_{\lambda'}\tilde\psi_{\lambda'}
    */
    void decompose_t(const InfiniteVector<double, Index>& c, const int jmin,
		     InfiniteVector<double, Index>& v) const;

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

  };
}

// include implementation
#include <interval/periodic.cpp>

#endif

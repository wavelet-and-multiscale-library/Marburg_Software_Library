// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PERIODIC_H
#define _WAVELETTL_PERIODIC_H

#include <iostream>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>

#include <interval/i_index.h>

using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    Template class for a periodic, biorthogonal wavelet basis on the unit interval [0,1],
    derived from a biorthogonal wavelet basis on R (which is specified as a
    template parameter).

    The periodized scaling functions (and analogously the wavelets) look like

      phi^per_{j,k}(x) = \sum_{l\in\mathbb Z} phi_{j,k}(x+l)

    References:
    [D] Daubechies,
        Ten Lectures On Wavelets, pp. 304ff
  */
  template <class RBASIS>
  class PeriodicBasis
  {
  public:
    //! default constructor
    PeriodicBasis();
    
    //! coarsest possible level j0
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef IntervalIndex<PeriodicBasis<RBASIS> > Index;

    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

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
      critical Sobolev regularity for the primal generators/wavelets
    */
    static inline double primal_regularity() { return RBASIS::primal_regularity(); }

    /*!
      degree of polynomial reproduction for the primal generators/wavelets
    */
    static unsigned int primal_polynomial_degree() { return RBASIS::primal_polynomial_degree(); }

    /*!
      number of vanishing moments for the primal wavelets
    */
    static inline unsigned int primal_vanishing_moments() { return RBASIS::primal_vanishing_moments(); }

    //! bounds for the generator indices
    inline const int DeltaLmin() const { return 0; }
    inline const int DeltaRmax(const int j) const { return (1<<j) - 1; }

    //! bounds for the wavelet indices
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j) - 1; }
    
    //! size of Delta_j
    inline const int Deltasize(const int j) const { return 1<<j; }
    
    //! size of Nabla_j
    inline const int Nablasize(const int j) const { return 1<<j; }

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
    //! an instance of the corresponding basis on R
    RBASIS rbasis;

    //! coarsest possible level
    int j0_;
  };
}

// include implementation
#include <interval/periodic.cpp>

#endif

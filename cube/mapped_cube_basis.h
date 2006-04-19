// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_MAPPED_CUBE_BASIS_H
#define _WAVELETTL_MAPPED_CUBE_BASIS_H

#include <geometry/chart.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>

using MathTL::Chart;

namespace WaveletTL
{
  /*!
    Template class for wavelet bases on the image of the d-dimensional (hyper)cube
    under a smooth parametrization kappa:\Box:=[0,1]^d -> R^m.
    The resulting wavelets are defined via

    \psi(x) = psi^\Box(.)/sqrt(g(.)) (kappa^{-1}(x)),

    where g() is the square root of the corresponding Gramian determinant.
    The reason for this normalization is that we require

    ||\psi||_{L_2(kappa(\Box))} = ||\psi^\Box||_{L_2(\Box)}.

    For each of the 2*d facets, you can specify the orders s_i, sT_i
    of the Dirichlet boundary conditions of the primal and dual basis
    in the normal direction (this is tailored for the DSBasis constructor...).
  */
  template <class IBASIS, unsigned int DIM_d = 2, unsigned int DIM_m = 2>
  class MappedCubeBasis
    : public CubeBasis<IBASIS,DIM_d>
  {
  public:
    /*!
      default constructor (no b.c.'s, identity chart),
      this is essentially used for testing
    */
    MappedCubeBasis();

    /*!
      constructor with a chart and specified boundary condition orders
      i-th direction at x=0 <-> index 2*i
      i-th direction at x=1 <-> index 2*i+1
    */
    MappedCubeBasis(const Chart<DIM_d,DIM_m>* kappa,
		    const FixedArray1D<int,2*DIM_d>& s,
		    const FixedArray1D<int,2*DIM_d>& sT);

    /*!
      constructor with a chart and specified boundary condition orders
      i-th direction at x=0 <-> index 2*i
      i-th direction at x=1 <-> index 2*i+1
    */
    MappedCubeBasis(const Chart<DIM_d,DIM_m>* kappa,
		    const FixedArray1D<int,2*DIM_d>& s);


    //! destructor
    ~MappedCubeBasis();

    //! wavelet index class
    typedef CubeIndex<IBASIS,DIM_d,MappedCubeBasis<IBASIS,DIM_d,DIM_m> > Index;

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
    //! the parametrization kappa: [0,1]^d -> R^m
    const Chart<DIM_d,DIM_m>* kappa_;

    //! flag for deletion of kappa
    bool delete_kappa;
  };
}

#include <cube/mapped_cube_basis.cpp>

#endif

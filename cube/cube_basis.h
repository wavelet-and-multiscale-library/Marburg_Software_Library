// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_BASIS_H
#define _WAVELETTL_CUBE_BASIS_H

#include <list>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>
#include <cube/cube_index.h>

// for convenience, include also some functionality
#include <cube/cube_support.h>

using std::list;
using MathTL::FixedArray1D;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    Template class for tensor product wavelet bases on the
    d-dimensional (hyper)cube [0,1]^d.
    For each of the 2*d facets, you can specify the orders s_i, sT_i
    of the Dirichlet boundary conditions of the primal and dual basis
    in the normal direction (this is tailored for the DSBasis constructor...).
  */
  template <class IBASIS, unsigned int DIM = 2>
  class CubeBasis
  {
  public:
    //! default constructor (no b.c.'s)
    CubeBasis();

    /*!
      constructor with specified boundary condition orders
      i-th direction at x=0 <-> index 2*i
      i-th direction at x=1 <-> index 2*i+1
    */
    CubeBasis(const FixedArray1D<int,2*DIM>& s,
	      const FixedArray1D<int,2*DIM>& sT);

    /*!
      constructor with specified Dirichlet boundary conditions for
      the primal functions, the dual functions will be constructed to
      fulfill free b.c.'s
    */
    CubeBasis(const FixedArray1D<bool,2*DIM>& bc);

    //! destructor
    ~CubeBasis();

    //! interval basis
    typedef IBASIS IntervalBasis;

    //! coarsest possible level j0
    inline const int j0() const { return j0_; }

    //! wavelet index class
    typedef CubeIndex<IBASIS,DIM,CubeBasis<IBASIS,DIM> > Index;

    //! read access to the bases
    const FixedArray1D<IBASIS*,DIM> bases() const { return bases_; }

    /*!
      geometric type of the support sets
    */
    typedef struct {
      int j;
      int a[DIM];
      int b[DIM];
    } Support;

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
    //! coarsest possible level j0
    int j0_;

    /*!
      the instances of the 1D bases (in general, we will of course
      need strictly less than DIM instances)
    */
    list<IBASIS*> bases_infact;

    //! for faster access, all relevant pointers to the 1D bases
    FixedArray1D<IBASIS*,DIM> bases_;
  };
}

#include <cube/cube_basis.cpp>

#endif

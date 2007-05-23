// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
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
#include <cube/cube_evaluate.h>

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
      constructor with specified boundary condition orders
      i-th direction at x=0 <-> index 2*i
      i-th direction at x=1 <-> index 2*i+1
    */
    CubeBasis(const FixedArray1D<int,2*DIM>& s);

    /*!
      constructor with specified Dirichlet boundary conditions for
      the primal functions, the dual functions will be constructed to
      fulfill free b.c.'s
    */
    CubeBasis(const FixedArray1D<bool,2*DIM>& bc);

    /*!
      constructor with precomputed instances of the 1D bases;
      in this case, the pointers are reused and
      not deleted at destruction time.
    */
    CubeBasis(const FixedArray1D<IBASIS*,DIM> bases);
    
    //! destructor
    ~CubeBasis();
    
    //! interval basis
    typedef IBASIS IntervalBasis;

    //! coarsest possible level j0
    inline const int j0() const { return j0_; }

    void set_jmax(const int jmax) {
      jmax_ = jmax;
      setup_full_collection();
    }

    //! wavelet index class
    typedef CubeIndex<IBASIS,DIM,CubeBasis<IBASIS,DIM> > Index;

    //! read access to the bases
    const FixedArray1D<IBASIS*,DIM>& bases() const { return bases_; }

    /*!
      geometric type of the support sets
    */
    typedef struct {
      int j;
      int a[DIM];
      int b[DIM];
    } Support;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return IBASIS::primal_regularity(); }
    
    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return IBASIS::primal_polynomial_degree(); }

    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }

    //! index of first generator on level j >= j0
    Index first_generator(const int j) const;
      
    //! index of last generator on level j >= j0
    Index last_generator(const int j) const;
      
    //! index of first wavelet on level j >= j0
    Index first_wavelet(const int j) const;
      
    //! index of last wavelet on level j >= j0
    Index last_wavelet(const int j) const;

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

    /*!
      For a given function, compute all integrals w.r.t. the primal
      or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
      - When integrating against the primal functions, the integrand has to be smooth
        to be accurately reproduced by the dual basis.
      - When integration against dual functions is specified,
        we integrate against the primal ones instead and multiply the resulting
        coefficients with the inverse of the primal gramian.

      Maybe a thresholding of the returned coefficients is helpful (e.g. for
      expansions of spline functions).
    */
    void
    expand
    (const Function<DIM>* f,
     const bool primal,
     const int jmax,
     InfiniteVector<double,Index>& coeffs) const;

    /*!
      helper function, integrate a smooth function f against a
      primal generator or wavelet
    */
    double
    integrate
    (const Function<DIM>* f,
     const Index& lambda) const;

    //! setup full collection of wavelets between j0_ and jmax_ as long as a jmax_ has been specified
    void setup_full_collection();

    //! collection of all wavelets between coarsest and finest level
    Array1D<Index> full_collection;

    //! number of wavelets between coarsest and finest level
    const int degrees_of_freedom() const { return full_collection.size(); };

    //! get the wavelet index corresponding to a specified number
    const inline Index* get_wavelet (const int number) const {
      return &full_collection[number];
    }

  protected:
    //! coarsest possible level j0
    int j0_;

    //! finest possible level j0
    int jmax_;

    /*!
      the instances of the 1D bases (in general, we will of course
      need strictly less than DIM instances)
    */
    list<IBASIS*> bases_infact;

    //! for faster access, all relevant pointers to the 1D bases
    FixedArray1D<IBASIS*,DIM> bases_;

    //! flag for deletion of the pointers
    bool delete_pointers;
  };
}

#include <cube/cube_basis.cpp>

#endif

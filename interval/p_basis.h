// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_P_BASIS_H
#define _WAVELETTL_P_BASIS_H

#include <algebra/matrix.h>
#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>
#include <interval/i_index.h>

using MathTL::Matrix;

namespace WaveletTL
{
  /*!
    Template class for the wavelet bases on the interval as introduced in [P].

    The primal generators are exactly those B-splines associated with the
    Schoenberg knot sequence

      t^j_{-d+1} = ... = t_0 = 0          (knot with multiplicity d at x=0)
      t^j_k = k * 2^{-j}, 1 <= k <= 2^j-1
      t^j_{2^j} = ... = t_{2^j+d-1} = 1   (knot with multiplicity d at x=1)

    i.e.

      B_{j,k}(x) = (t^j_{k+d}-t^j_k)[t^j_k,...,t^j_{k+d}](t-x)^{d-1}_+

    with supp(B_{j,k}(x) = [t^j_k, t^j_{k+d}].
    In other words, we have exactly

     d-1     left boundary splines  (k=-d+1,...,-1),
     2^j-d+1 inner splines          (k=0,...,2^j-d),
     d-1     right boundary splines (k=2^j-d+1,...,2^j-1)

    Since the primal CDF generators are centered around floor(d/2)=-ell_1,
    we perform an index shift by ell_1, i.e., we use the generators
 
      phi_{j,k}(x) = 2^{j/2} B_{j,k-ell_1}

    So, if no boundary conditions are imposed, the index of the leftmost generator
    will be 1-d+floor(d/2). See default constructor for possible b.c.'s
    
    References:
    [P] Primbs:
        Stabile biorthogonale Wavelet-Basen auf dem Intervall
	Dissertation, Univ. Duisburg-Essen, 2006
  */
  template <int d, int dT>
  class PBasis
  {
  public:
    /*!
      constructor
      
      At the moment, you may (only) specify the order of the primal (s) boundary conditions
      at the left and right end of the interval [0,1].
      For technical reasons (easing the setup of the dual generators),
      only the values si >= d-2 are currently implemented.
      In other words, linear spline wavelets may or may not fulfill Dirichlet b.c.'s.
      Higher order wavelet bases will always fulfill b.c.'s of order at least d-2.
      The dual wavelet basis will have no b.c.'s in either case and can reproduce the
      full range of polynomials of order dT.
    */
    PBasis(const int s0 = 0, const int s1 = 0);

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! freezing parameters
    inline const int ellT_l() const { return -ell1T<d,dT>()+s0+2-d; }
    inline const int ellT_r() const { return -ell1T<d,dT>()+s1+2-d; }
    
    //! wavelet index class
    typedef IntervalIndex<PBasis<d,dT> > Index;
    
    //! extremal generator indices
    inline const int DeltaLmin() const { return 1-d-ell1<d>()+s0; }
    inline const int DeltaRmax(const int j) const { return (1<<j)-1-ell1<d>()-s1; }

    //! boundary indices in \nabla_j
    inline const int Nablamin() const { return 0; }
    inline const int Nablamax(const int j) const { return (1<<j)-1; }
    
  protected:
    //! coarsest possible level
    int j0_;

    //! boundary condition orders at 0 and 1
    int s0, s1;

    //! general setup routine which is shared by the different constructors
    void setup();

    //! one instance of a CDF basis (for faster access to the primal and dual masks)
    CDFBasis<d,dT> cdf;

    //! single CDF moments \alpha_{m,r} := \int_{\mathbb R} x^r\phi(x-m)\,dx
    const double alpha(const int m, const unsigned int r) const;

     //! refinement coeffients of left dual boundary generators
    const double betaL(const int m, const unsigned int r) const;

   //! boundary blocks in Mj0
    Matrix<double> ML_, MR_;

    //! boundary blocks in Mj0T
    Matrix<double> MLT_, MRT_;

    //! Gramian matrices for the left and right generators (primal against unbiorth. dual)
    Matrix<double> GammaL, GammaR;
  };
}

#include <interval/p_basis.cpp>

#endif

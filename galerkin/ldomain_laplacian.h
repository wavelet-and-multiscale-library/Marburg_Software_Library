// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_LAPLACIAN_H
#define _WAVELETTL_LDOMAIN_LAPLACIAN_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>

#include <interval/spline_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor> class SplineBasis;
  template <int d, int dT> class SplineBasis<d,dT,DS_construction>;
 
  template <class IBASIS> class LDomainBasis;
  template <int d, int dT> class LDomainBasis<SplineBasis<d,dT,DS_construction> >;

  /*!
    This class models the stiffness matrix of the Laplacian
    with respect to a composite wavelet basis
    on the L--shaped domain in R^2.

    The general template class is incomplete!
  */
  template <class IBASIS>
  class LDomainLaplacian
    : public FullyDiagonalEnergyNormPreconditioner<typename LDomainBasis<IBASIS>::Index>
  {
  public:
    /*!
      make template argument accessible
    */
    typedef LDomainBasis<IBASIS> WaveletBasis;

    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

  };
  
  // template specialization for the case IBASIS==SplineBasis<d,dT,DS_construction>
  template <int d, int dT> class LDomainLaplacian<SplineBasis<d,dT,DS_construction> >;

  template <int d, int dT>
  class LDomainLaplacian<SplineBasis<d,dT,DS_construction> >
    : public FullyDiagonalEnergyNormPreconditioner<typename LDomainBasis<SplineBasis<d,dT,DS_construction> >::Index>
  {
  public:
    /*!
      the 1D wavelet basis class
    */
    typedef SplineBasis<d,dT,DS_construction> IntervalBasis;
    
    /*!
      the wavelet basis class
    */
    typedef LDomainBasis<IntervalBasis> WaveletBasis;
    
    /*!
      wavelet index class
    */
    typedef typename WaveletBasis::Index Index;

    /*!
      constructor from a given wavelet basis and a given right-hand side y
    */
    LDomainLaplacian(const WaveletBasis& basis,
		     const InfiniteVector<double,Index>& y);
    
    /*!
      set right-hand side y
    */
    void set_rhs(const InfiniteVector<double,Index>& y) const {
      y_ = y;
    }

    /*!
      evaluate the diagonal preconditioner D
    */
    double D(const Index& lambda) const { return sqrt(a(lambda,lambda)); }
    
    /*!
      evaluate the (unpreconditioned) bilinear form a
    */
    double a(const Index& lambda,
 	     const Index& nu) const;

    /*!
      evaluate the (unpreconditioned) right-hand side f
    */
    double f(const Index& lambda) const {
      return y_.get_coefficient(lambda);
    }

  protected:
    const WaveletBasis& basis_;
    
    // rhs, mutable to have 'const' method
    mutable InfiniteVector<double,Index> y_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };
}

#include <galerkin/ldomain_laplacian.cpp>

#endif

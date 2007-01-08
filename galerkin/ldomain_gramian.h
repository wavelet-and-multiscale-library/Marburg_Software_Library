// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_GRAMIAN_H
#define _WAVELETTL_LDOMAIN_GRAMIAN_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>

#include <interval/spline_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <galerkin/ldomain_equation.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;

namespace WaveletTL
{
  template <int d, int dT, SplineBasisFlavor flavor> class SplineBasis;
  template <int d, int dT> class SplineBasis<d,dT,DS_construction>;
 
  template <class IBASIS> class LDomainBasis;
  template <int d, int dT> class LDomainBasis<SplineBasis<d,dT,DS_construction> >;

  /*!
    This class models the Gramian matrix of a composite wavelet basis
    on the L--shaped domain in R^2.
  */
  template <class IBASIS>
  class LDomainGramian
    : public LDomainEquation<IBASIS>
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

    /*!
      constructor from an identity bvp
    */
    LDomainGramian(IdentityBVP<2>* bvp)
      : LDomainEquation<IBASIS>(bvp, false)
    {
    }

    /*!
      (half) order t of the operator
    */
    double operator_order() const { return 0.; }
    
    /*!
      evaluate the diagonal preconditioner D (we don't have any)
    */
    double D(const Index& lambda) const { return 1; }
  };
  
  // template specialization for the case IBASIS==SplineBasis<d,dT,DS_construction>
  template <int d, int dT>
  class LDomainGramian<SplineBasis<d,dT,DS_construction> >
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
    LDomainGramian(const WaveletBasis& basis,
 		   const InfiniteVector<double,Index>& y);
    
  protected:
    const WaveletBasis& basis_;
    
    // rhs, mutable to have 'const' method
    mutable InfiniteVector<double,Index> y_;
    
    // estimates for ||A|| and ||A^{-1}||
    mutable double normA, normAinv;
  };
}

#include <galerkin/ldomain_gramian.cpp>

#endif

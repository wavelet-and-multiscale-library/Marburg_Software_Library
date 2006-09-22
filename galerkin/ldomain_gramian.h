// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_GRAMIAN_H
#define _WAVELETTL_LDOMAIN_GRAMIAN_H

#include <set>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <numerics/bvp.h>

#include <galerkin/galerkin_utils.h>
#include <galerkin/infinite_preconditioner.h>
#include <galerkin/ldomain_equation.h>

using MathTL::FixedArray1D;
using MathTL::EllipticBVP;

namespace WaveletTL
{
  template <class IBASIS> class LDomainBasis;

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
}

#include <galerkin/ldomain_gramian.cpp>

#endif

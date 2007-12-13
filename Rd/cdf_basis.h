// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_BASIS_H
#define _WAVELETTL_CDF_BASIS_H

#include <Rd/cdf_mask.h>
#include <Rd/r_basis.h>

namespace WaveletTL
{
  /*!
    Biorthogonal wavelet basis on R as introduced in [CDF].
    In addition to the generic routines from RBasis,
    point evaluation of the primal generators and wavelets is provided.

    References:
    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
   */
  template <int d, int dt>
  class CDFBasis
    : public virtual RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >
  {
  public:
    // use generic point evaluation routines on dyadic grids from RBasis
    using RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::evaluate;

    // provide masks
    using RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::primal_mask;
    using RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::dual_mask;

    //! point evaluation of (derivatives of) primal generators and wavelets
    double evaluate(const unsigned int derivative,
		    const RIndex& lambda,
		    const double x) const;
  };
}

#include <Rd/cdf_basis.cpp>

#endif

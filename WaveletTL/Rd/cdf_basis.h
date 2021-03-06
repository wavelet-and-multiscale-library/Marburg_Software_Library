// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_BASIS_H
#define _WAVELETTL_CDF_BASIS_H

#include <utils/array1d.h>
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
    // provide masks
    using RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::primal_mask;
    using RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::dual_mask;

    // use generic point evaluation routines on dyadic grids from RBasis
    using RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::evaluate;

    // use support routines from RBasis
    using RBasis<CDFRefinementMask_primal<d>, CDFRefinementMask_dual<d, dt> >::support;

    //! point evaluation of (derivatives of) primal generators and wavelets
    static double evaluate(const unsigned int derivative,
			   const RIndex& lambda,
			   const double x);

    /*!
      point evaluation of (derivatives) of a single primal generator
      or wavelet \psi_\lambda at several points simultaneously
    */
    static void evaluate
    (const unsigned int derivative,
     const RIndex& lambda,
     const Array1D<double>& points, Array1D<double>& values);
    
    //integration in [-\infty,\infty] of the generator or wavelet \psi_\lambda 
    static double integrate(const RIndex& lambda);
  };
}

#include <Rd/cdf_basis.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_QUARKLET_FRAME_H
#define _WAVELETTL_QUARKLET_FRAME_H

#include <utils/array1d.h>
#include <Rd/quarklet_mask.h>
#include <Rd/r_frame.h>

namespace WaveletTL
{
  /*!
    Quarklet Frame on R based on CDF-Wavelets (see [CDF]).
    In addition to the generic routines from RFrame,
    point evaluation of the primal generators and quarklets is provided.

    References:
    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
   */
  template <int d, int dt>
  class QuarkletFrame
    : public virtual RFrame<QuarkletRefinementMask_primal<d>, QuarkletRefinementMask_dual<d, dt> >
  {
  public:
    // provide masks
    using RFrame<QuarkletRefinementMask_primal<d>, QuarkletRefinementMask_dual<d, dt> >::primal_mask;
    using RFrame<QuarkletRefinementMask_primal<d>, QuarkletRefinementMask_dual<d, dt> >::dual_mask;

    // use generic point evaluation routines on dyadic grids from RFrame
    using RFrame<QuarkletRefinementMask_primal<d>, QuarkletRefinementMask_dual<d, dt> >::evaluate;

    // use support routines from RBasis
    using RFrame<QuarkletRefinementMask_primal<d>, QuarkletRefinementMask_dual<d, dt> >::support;
    
    using RFrame<QuarkletRefinementMask_primal<d>, QuarkletRefinementMask_dual<d, dt> >::reconstruct_1;

    //! point evaluation of (derivatives of) primal generators and quarklets
    static double evaluate(const unsigned int derivative,
			   const RQIndex& lambda,
			   const double x);

    /*!
      point evaluation of (derivatives) of a single primal generator
      or quarklet \psi_\lambda at several points simultaneously
    */
    static void evaluate
    (const unsigned int derivative,
     const RQIndex& lambda,
     const Array1D<double>& points, Array1D<double>& values);
    
    //integration in [-\infty,\infty] of the quark or quarklet \psi_\lambda
    //easy to expand to derivatives (not necessary a.t.m.)
    static double integrate(const RQIndex& lambda);
    
    private:
      //Instance of QuarkletFrame
      //QuarkletFrame<d,dt> q_frame;
  };
  
}

#include <Rd/quarklet_frame.cpp>

#endif

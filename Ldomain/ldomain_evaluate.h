// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_EVALUATE_H
#define _WAVELETTL_LDOMAIN_EVALUATE_H

#include <utils/array1d.h>
#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <Ldomain/ldomain_basis.h>
#include <interval/spline_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <class IBASIS> class LDomainBasis;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of the L-shaped domain

    We redirect this to a member function of basis.
  */
  template <class IBASIS>
  Array1D<SampledMapping<2> >
  evaluate(const LDomainBasis<IBASIS>& basis,
	   const typename LDomainBasis<IBASIS>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    return basis.evaluate(lambda, primal, resolution);
  }
   
  /*!
    Evaluate an arbitrary linear combination of primal/dual wavelets
    on a dyadic subgrid of the L-shaped domain

    We redirect this to a member function of basis.
  */
  template <class IBASIS>
  Array1D<SampledMapping<2> >
  evaluate(const LDomainBasis<IBASIS>& basis,
	   const InfiniteVector<double, typename LDomainBasis<IBASIS>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    return basis.evaluate(coeffs, primal, resolution);
  }
}

#endif

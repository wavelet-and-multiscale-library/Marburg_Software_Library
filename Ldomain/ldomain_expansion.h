// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_EXPANSION_H
#define _WAVELETTL_LDOMAIN_EXPANSION_H

#include <utils/function.h>
#include <algebra/infinite_vector.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <class IBASIS> class LDomainBasis;

  /*!
    helper function, integrate a smooth function f against a
    primal Ldomain generator or wavelet
  */
  template <class IBASIS>
  double integrate(const Function<2>* f,
		   const LDomainBasis<IBASIS>& basis,
		   const typename LDomainBasis<IBASIS>::Index& lambda);

  /*!
    For a given function, compute all integrals w.r.t. the primal
    generators/wavelets (primal=true)
    or approximate the expansion coefficients w.r.t. the primal basis
    (primal=false)

    remark: at the moment, "primal" is a dummy variable yet.
  */
  template <class IBASIS>
  void
  expand(const Function<2>* f,
	 const LDomainBasis<IBASIS>& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double, typename LDomainBasis<IBASIS>::Index>& coeffs);
}

#include <Ldomain/ldomain_expansion.cpp>

#endif

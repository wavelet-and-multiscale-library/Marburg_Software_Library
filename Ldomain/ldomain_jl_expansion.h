// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_JL_EXPANSION_H
#define _WAVELETTL_LDOMAIN_JL_EXPANSION_H

#include <utils/function.h>
#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <Ldomain/ldomain_jl_basis.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;
using MathTL::Function;

namespace WaveletTL
{
  typedef LDomainJLBasis::Index Index;

  /*!
    helper function, integrate a smooth function f against a
    primal Ldomain generator or wavelet
  */
  double integrate(const Function<2>* f,
		   const LDomainJLBasis& basis,
		   const Index& lambda);

  /*!
    For a given function, compute all integrals w.r.t. the primal
    generators/wavelets (primal=true)
    or approximate the expansion coefficients w.r.t. the primal basis
    (primal=false)
  */
  void
  expand(const Function<2>* f,
	 const LDomainJLBasis& basis,
	 const bool primal,
	 const int jmax,
	 InfiniteVector<double,Index>& coeffs);
}

#include <Ldomain/ldomain_jl_expansion.cpp>

#endif

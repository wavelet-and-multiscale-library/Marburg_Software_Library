// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_BOUNDARY_GRAMIAN_H
#define _WAVELETTL_BOUNDARY_GRAMIAN_H

#include <algebra/matrix.h>
#include <Rd/cdf_utils.h>
#include <Rd/cdf_basis.h>

using MathTL::Matrix;

namespace WaveletTL
{
  /*!
    Assume that phi, phiT are refinable functions on R,
    with compact supports [l1,l2]\subset [l1T,l2T].
    It is known that the inner products <phi(.-k),phiT(.-l)> over R can be
    computed exactly, since they correspond to the point values at the integers
    of a special compactly supported refinable function, see the
    template class DMMask for their refinement coefficients.

    The following routine computes the gramian matrix
      Gamma^L = <phi^L_{0,k},phiT^L_{0,l}>_{k,l=k0,...,ell-1}
    between boundary functions that satisfy a refinement relation of the form
      (phi^L_{0,k})_{k=k0,...,ell-1}=M_L^T (phi^L_{0,k}(2.))_{k=k0,...}

  */
  template <class MASK1, class MASK2>
  void
  compute_boundary_gramian(Matrix<double>& Gamma);
}

#include <interval/boundary_gramian.cpp>

#endif

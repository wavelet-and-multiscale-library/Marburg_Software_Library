// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_SPLINE_BASIS_H
#define _WAVELETTL_LDOMAIN_SPLINE_BASIS_H

#include <interval/spline_basis.h>

namespace WaveletTL
{
  /*!
    Template class for composite wavelet bases over the L-shaped domain
      (-1,1)^2 \ (0,1)^2
    with complementary (!) b.c.'s at the outer domain boundary.
    The primal generators are essentially glued B-splines from the
    local subpatches.
    
    References:
    [DS] Dahmen, Schneider:
         Composite Wavelet Bases for Operator Equations
	 Math. Comput. 68 (1999), 1533-1567
  */
  template <int d, int dT>
  class LDomainSplineBasis
  {
  public:
    //! default constructor
    LDomainSplineBasis();

  protected:
    // internal spline bases
    SplineBasis<d,dT> basis1D_01, basis1D_10;
  };
}

#include <Ldomain/ldomain_spline_basis.cpp>

#endif

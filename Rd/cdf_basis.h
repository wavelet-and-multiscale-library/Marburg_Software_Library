// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_BASIS_H
#define _WAVELETTL_CDF_BASIS_H

#include <Rd/cdf_mask.h>
#include <Rd/cdf_utils.h>
#include <Rd/r_basis.h>

namespace WaveletTL
{
  // since C++ does not allow typedef template, we use this simple wrapper class:
  template <int d, int dt>
  class CDFBasis
    : public RBasis<CDFMask_primal<d>, CDFMask_dual<d, dt> >
  {
  };
}

#endif

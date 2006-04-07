// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_BASIS_H
#define _WAVELETTL_LDOMAIN_BASIS_H

#include <list>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/fixed_array1d.h>
#include <utils/multiindex.h>
#include <cube/cube_evaluate.h>

using std::list;
using MathTL::FixedArray1D;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  /*!
    Template class for composite wavelet bases over the L-shaped domain
      [-1,1]^2 \ (0,1]^2
    with homogeneous b.c.'s for the primal basis and free b.c.'s for the
    dual basis.
  */
  template <class IBASIS>
  class LDomainBasis
  {
  public:
  };
}

#include <Ldomain/ldomain_basis.cpp>

#endif

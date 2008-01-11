// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_RING_BASIS_H
#define _WAVELETTL_RING_BASIS_H

#include <iostream>

#include <algebra/vector.h>
#include <algebra/infinite_vector.h>
#include <utils/array1d.h>
#include <utils/function.h>
#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <interval/spline_basis.h>

using MathTL::Vector;
using MathTL::InfiniteVector;
using MathTL::Array1D;
using MathTL::Function;

namespace WaveletTL
{
  /*!
    Template class for a wavelet basis on the ring-shaped domain
      R = {(x,y) : r_0 <= ||(x,y)|| <= r_1 }
    with complementary boundary conditions of order s0 at the interior boundary
    and of order s1 at the outer boundary.
  */
  template <int d, int dt, int s0=0, int s1=0>
  class RingBasis
  {
  public:
    /*!
      coarsest possible level j0; we may assume that the interval basis in radial direction
      has a larger minimal level than the corresponding periodic basis
    */
    static const int j0() {
      return SplineBasis<d,dt,P_construction,s0,s1,0,0>::j0();
    }

  protected:
    //! an instance of the periodic basis (angular direction)
    PeriodicBasis<CDFBasis<d,dt> > periodic_basis;

    //! an instance of the nonperiodic basis (radial direction)
    SplineBasis<d,dt,P_construction,s0,s1,0,0> nonperiodic_basis;
  };
}

// include implementation
#include <ring/ring_basis.cpp>

#endif

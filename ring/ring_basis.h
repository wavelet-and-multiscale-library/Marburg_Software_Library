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
#include <geometry/sampled_mapping.h>

#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <interval/spline_basis.h>
#include <ring/ring_index.h>

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
    //! default constructor, also sets r0 and r1
    RingBasis(const double r0 = 1.0, const double r1 = 2.0);
    
    /*!
      coarsest possible level j0; we may assume that the interval basis in radial direction
      has a larger minimal level than the corresponding periodic basis
    */
    static const int j0() {
      return SplineBasis<d,dt,P_construction,s0,s1,0,0>::j0();
    }

    //! wavelet index class
    typedef RingIndex<d,dt,s0,s1> Index;

    //! index of first generator on level j >= j0
    static Index first_generator(const int j);

    //! index of last generator on level j >= j0
    static Index last_generator(const int j);

    //! index of first wavelet on level j >= j0
    static Index first_wavelet(const int j);

    //! index of last wavelet on level j >= j0
    static Index last_wavelet(const int j);


    /*!
      Evaluate a single primal/dual generator or wavelet \psi_\lambda
      on a "dyadic" subgrid of the ring, pulled back to the unit square
    */
    SampledMapping<2> evaluate(const Index& lambda,
			       const int resolution) const;
    
    /*!
      Evaluate an arbitrary linear combination of primal/dual wavelets
    */
    SampledMapping<2> evaluate(const InfiniteVector<double,Index>& coeffs,
 			       const int resolution) const;
    

    /*!
      For a given function, compute all integrals w.r.t. the primal
      or dual generators/wavelets \psi_\lambda with |\lambda|\le jmax.
    */
    void
    expand
    (const Function<2>* f,
     const bool primal,
     const int jmax,
     InfiniteVector<double, Index>& coeffs) const;

    /*!
      helper function, integrate a smooth function f against a
      primal generator or wavelet
    */
    double
    integrate
    (const Function<2>* f,
     const Index& lambda) const;

    /*!
      helper function, integrate two primal generators or wavelets
      against each other (for the Gramian)
    */
    double
    integrate
    (const Index& lambda,
     const Index& mu) const;

  protected:
    //! an instance of the periodic basis (angular direction)
    PeriodicBasis<CDFBasis<d,dt> > basis0;

    //! an instance of the nonperiodic basis (radial direction)
    SplineBasis<d,dt,P_construction,s0,s1,0,0> basis1;

    //! radii
    double r0_, r1_;
  };
}

// include implementation
#include <ring/ring_basis.cpp>

#endif

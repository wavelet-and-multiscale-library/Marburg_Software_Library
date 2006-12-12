// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_SPLINE_BASIS_H
#define _WAVELETTL_LDOMAIN_SPLINE_BASIS_H

#include <algebra/vector.h>
#include <interval/spline_basis.h>

using namespace MathTL;

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

    //! coarsest possible level
    inline const int j0() const { return j0_; }

    //! size_type, for convenience
    typedef Vector<double>::size_type size_type;

    //! interval basis
    typedef SplineBasis<d,dT> IntervalBasis;
    
    //! geometric type of the support sets
    typedef struct {
      int j;       // granularity
      int xmin[3];
      int xmax[3];
      int ymin[3];
      int ymax[3];
    } Support;
    
    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return IntervalBasis::primal_regularity(); } // dirty, we should use max(1.5,~) instead
    
    //! degree of polynomial reproduction for the primal generators/wavelets
    static unsigned int primal_polynomial_degree() { return IntervalBasis::primal_polynomial_degree(); }
    
    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return IntervalBasis::primal_vanishing_moments(); }
    
    

  protected:
    // internal spline bases
    SplineBasis<d,dT> basis1D_01, basis1D_10, basis1D_11;
    
    // coarsest level
    int j0_;
  };
}

#include <Ldomain/ldomain_spline_basis.cpp>

#endif

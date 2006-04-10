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

#include <Ldomain/ldomain_index.h>

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
    //! default constructor
    LDomainBasis();

    //! interval basis
    typedef IBASIS IntervalBasis;

    //! wavelet index class
    typedef LDomainIndex<IBASIS> Index;

    //! critical Sobolev regularity for the primal generators/wavelets
    static double primal_regularity() { return IBASIS::primal_regularity(); }
    
    //! number of vanishing moments for the primal wavelets
    static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }

    //! number of vanishing moments for the dual wavelets
    static unsigned int dual_vanishing_moments() { return IBASIS::dual_vanishing_moments(); }

  protected:
    //! the two interval wavelet bases involved
    IntervalBasis basis01, basis10;
  };
}

#include <Ldomain/ldomain_basis.cpp>

#endif

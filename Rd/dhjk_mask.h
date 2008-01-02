// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DHJK_MASK_H
#define _WAVELETTL_DHJK_MASK_H

#include <iostream>
#include <algebra/infinite_vector.h>
#include <algebra/fixed_matrix.h>

using MathTL::InfiniteVector;
using MathTL::FixedMatrix;

namespace WaveletTL
{
  /*!
    primal mask of the biorthogonal refinable multigenerator constructed in

    [DHJK] Dahmen, Han, Jia, Kunoth: Biorthogonal Multiwavelets on the Interval:
           Cubic Hermite Splines.  Constr. Approx. (2000) 16:221-259
  */
  class DHJKMask_primal
    : public InfiniteVector<FixedMatrix<double, 2> >
  {
  public:
    /*!
      Default constructor
      constructs the refinement mask
    */
    DHJKMask_primal();

    /*!
      number of compononents of the multi-wavelet basis
    */
    static const unsigned int number_of_components = 2;

    /*!
      Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    */
    static unsigned int Strang_Fix_order() { return 4; }

    /*!
      critical Sobolev regularity
    */
    static double regularity() { return 1.5; }
  };

  /*!
    dual mask of the biorthogonal refinable multigenerator constructed in

    [DHJK] Dahmen, Han, Jia, Kunoth: Biorthogonal Multiwavelets on the Interval:
           Cubic Hermite Splines.  Constr. Approx. (2000) 16:221-259
  */
  class DHJKMask_dual
    : public InfiniteVector<FixedMatrix<double, 2> >
  {
  public:
    /*!
      Default constructor
      constructs the refinement mask
    */
    DHJKMask_dual();

    /*!
      number of compononents of the multi-wavelet basis
    */
    static const unsigned int number_of_components = 2;

    /*!
      Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    */
    static unsigned int Strang_Fix_order() { return 2; }

    /*!
      critical Sobolev regularity
      (estimate, cf. [DHJK, Prop. 3.2])
    */
    static double regularity() { return 0.824926; }
  };
}

// include implementation
#include <Rd/dhjk_mask.cpp>

#endif // _WAVELETTL_DHJK_MASK_H

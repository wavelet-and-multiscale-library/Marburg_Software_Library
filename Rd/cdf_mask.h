// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_MASK_H
#define _WAVELETTL_CDF_MASK_H

#include <iostream>
#include <algebra/multi_laurent_polynomial.h>

using MathTL::MultivariateLaurentPolynomial;

namespace WaveletTL
{
  /*!
    primal mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d>
  class CDFMask_primal
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    CDFMask_primal();

    /*!
      Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    */
    static unsigned int Strang_Fix_order() { return d; }

    /*!
      critical Sobolev regularity
    */
    static double regularity() { return d - 0.5; }
  };

  /*!
    dual mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d, int dT>
  class CDFMask_dual
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    CDFMask_dual();

    /*!
      Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    */
    static unsigned int Strang_Fix_order() { return dT; }

    /*!
      critical Sobolev regularity
      (at least a crude lower estimate, this routine will only cf. [CDF])

    */
    static double regularity() { return 1.0; }
  };
}

// include implementation
#include <Rd/cdf_mask.cpp>

#endif

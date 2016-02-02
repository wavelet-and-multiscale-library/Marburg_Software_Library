// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_QUARKLET_MASK_H
#define _WAVELETTL_QUARKLET_MASK_H

#include <iostream>
#include <algebra/multi_laurent_polynomial.h>
#include <Rd/r_q_mask.h>
#include <Rd/cdf_utils.h>

using MathTL::MultivariateLaurentPolynomial;

namespace WaveletTL
{
  /*!
    primal mask of the biorthogonal refinable functions constructed in [CDF]

    References:
    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d>
  class QuarkletMask_primal
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    QuarkletMask_primal();

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
    primal mask of the biorthogonal refinable functions constructed in [CDF]

    References:
    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d>
  class QuarkletRefinementMask_primal
    : public RQRefinementMask<d+1, -(d/2)> // ell_1 = -(d/2)
  {
  public:
    //! default constructor, sets up the mask coefficients
    QuarkletRefinementMask_primal();
    
    //! Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    static unsigned int Strang_Fix_order() { return d; }
    
    //! critical L_2 Sobolev regularity
    static double regularity() { return d - 0.5; }  
  };

  /*!
    dual mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d, int dt>
  class QuarkletMask_dual
    : public virtual MultivariateLaurentPolynomial<double, 1>
  {
  public:
    //! default constructor, sets up the mask coefficients
    QuarkletMask_dual();

    //! Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    static unsigned int Strang_Fix_order() { return dt; }

    //! critical L_2 Sobolev regularity (at least a crude lower estimate)
    static double regularity() { return 1.0; }
  };

  /*!
    dual mask of the biorthogonal refinable functions constructed in

    [CDF] Cohen, Daubechies, Feauveau: Biorthogonal bases of compactly supported wavelets,
          Comm. Pure Appl. Math. 45(1992), 485--560.
  */
  template <int d, int dt>
  class QuarkletRefinementMask_dual
    : public virtual RQRefinementMask<d+2*dt-1, -(d/2)-dt+1> // ellT1 = -(d/2)-dt+1
  {
  public:
    //! default constructor, sets up the mask coefficients
    QuarkletRefinementMask_dual();

    //! Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    static unsigned int Strang_Fix_order() { return dt; }

    //! critical L_2 Sobolev regularity (at least a crude lower estimate)
    static double regularity() { return 1.0; }
  };
}

// include implementation
#include <Rd/quarklet_mask.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DM_MASK_H
#define _WAVELETTL_DM_MASK_H

#include <algebra/laurent_polynomial.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    According to [DM], for refinable functions phi_0, phi_1, the function

      F(t) = \int_{\mathbb R}\phi_0(x)\phi_1(x-t)\,dx

    is again refinable. This template class computes the mask of F from those
    of phi_0 and phi_1. You can directly use the class as a template parameter
    of RefinableFunction<>.

    TODO: also handle the case of derivatives and the restriction to [0,1]

    References:
    [DM] W. Dahmen and C. Micchelli,
         Using the refinement equation for evaluating integrals of wavelets,
	 SIAM J. Numer. Anal. 30(1993), 507-537
  */
  template <class MASK0, class MASK1>
  class DMMask
    : public LaurentPolynomial<double>
  {
  public:
    //! default constructor, set up the mask
    DMMask();
  };
}

// include implementation
#include <Rd/dm_mask.cpp>

#endif

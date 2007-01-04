// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DM_MASK_H
#define _WAVELETTL_DM_MASK_H

#include <algebra/multi_laurent_polynomial.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    According to [DM], for refinable functions phi_0, phi_1, the function

      F(t) = \int_{\mathbb R}\phi_0(x)\phi_1(x-t)\,dx

    is again refinable. This template class computes the mask of F from those
    of phi_0 and phi_1. You can directly use the class as a template parameter
    of RefinableFunction<>.

    Make sure that the template parameters refer to univariate functions.

    References:
    [DM] W. Dahmen and C. Micchelli,
         Using the refinement equation for evaluating integrals of wavelets,
	 SIAM J. Numer. Anal. 30(1993), 507-537
  */
  template <class MASK0, class MASK1>
  class DMMask1
    : public MultivariateLaurentPolynomial<double, 1>
  {
  public:
    //! default constructor, set up the mask
    DMMask1();
  };

  /*!
    According to [DM], for refinable functions phi_0, phi_1, phi_2, the function

      F(t_1, t_2) = \int_{\mathbb R}\phi_0(x)\phi_1(x-t_1)\phi_2(x-t_2)\,dx

    is again refinable. This template class computes the mask of F from those
    of phi_0, phi_1 and phi_2. You can directly use the class as a template parameter
    of MultivariateRefinableFunction<2, *>.

    References:
    [DM] W. Dahmen and C. Micchelli,
         Using the refinement equation for evaluating integrals of wavelets,
	 SIAM J. Numer. Anal. 30(1993), 507-537
  */
  template <class MASK0, class MASK1, class MASK2>
  class DMMask2
    : public MultivariateLaurentPolynomial<double, 2>
  {
  public:
    //! default constructor, set up the mask
    DMMask2();
  };
}

// include implementation
#include <Rd/dm_mask.cpp>

#endif

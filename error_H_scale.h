// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_ERROR_H_SCALE_H
#define _FRAMETL_ERROR_H_SCALE_H

namespace FrameTL
{

  /*!
    \file error_H_scale.h
    Routines for the computation of the \f$L_2\f$- and \f$H^1\f$ errors of approximate
    solutions of the Poisson equation on the interval and in the L-shaped domain.
   */

  /*!\brief  This function computes the \f$L_2\f$-norm or \f$H^1\f$-seminorm of the difference between the frame
    expansion given by the coefficients in coeffs and a given function in the unit interval.
    
    The last parameter of this routine has to be the function, or, in case of the \f$H^1\f$-seminorm,
    it has to be the derivative of the function.
  */
  template <class IBASIS>
  double
  error_H_scale_interval (const int order,
			  const AggregatedFrame<IBASIS,1,1>& frame,
			  const InfiniteVector<double, typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
			  const Function<1>& f);

  /*!\brief This function computes the \f$L_2\f$-norm or \f$H^1\f$-seminorm of the difference between the frame
    expansion given by the coefficients in coeffs and a given function in the L-shaped domain
    \f$[-1,1]^2 \ [0,1)^2\f$.

    The last parameter of this routine has to be the function, or, in case of
    the \f$H^1\f$-seminorm, it has to be the gradient of the function.
   */
  template <class IBASIS>
  double
  error_H_scale_Lshaped (const int order,
			 const AggregatedFrame<IBASIS,2,2>& frame,
			 const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
 			 const Function<2>& f);
  /*! \brief Computes an approximation to the \f$H^1\f$-norm of a function in the L-shaped domain
    \f$[-1,1]^2 \ [0,1)^2\f$.

    The parameter gradient is assumed to be the gradient of the function.
  */
  double
  H_1_semi_norm_Lshaped(const Function<2>& gradient);

}

#include <error_H_scale.cpp>

#endif

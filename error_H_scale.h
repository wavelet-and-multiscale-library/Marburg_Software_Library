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
  template <class IBASIS>
  double
  error_H_scale_interval (const int order,
			  const AggregatedFrame<IBASIS,1,1>& frame,
			  const InfiniteVector<double, typename AggregatedFrame<IBASIS,1,1>::Index>& coeffs,
			  const Function<1>& f);

  template <class IBASIS>
  double
  error_H_scale_Lshaped (const int order,
			 const AggregatedFrame<IBASIS,2,2>& frame,
			 const InfiniteVector<double, typename AggregatedFrame<IBASIS,2,2>::Index>& coeffs,
 			 const Function<2>& f);
  /*!
    computes an approximation to the H^1 norm of a function on the L shaped domain
    [-1,1]^2 \ [0,1)^2,
    f is assumed to be the gradient of the function
  */
  double
  H_1_semi_norm_Lshaped(const Function<2>& f);

}

#include <error_H_scale.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_JL_UTILS_H
#define _WAVELETTL_JL_UTILS_H

#include <utils/array1d.h>

using MathTL::Array1D;

namespace WaveletTL
{
  // some standalone utility functions for JLBasis and variants thereof

  /*!
    point evaluation of (derivatives) of a single primal [JL] generator
    or wavelet \psi_\lambda;
    this routine works for _all_ Hermite spline wavelets on R
  */
  double evaluate(const unsigned int derivative,
		  const int j, const int e, const int c, const int k,
 		  const double x);
  
  /*!
    point evaluation of 0-th and first derivative of a single primal [JL] generator
    or wavelet \psi_\lambda at several points simultaneously;
    without a temporary JLBasis::Index object, the routine works for _all_ Hermite spline
    wavelets on R
  */
  void evaluate(const int j, const int e, const int c, const int k,
		const Array1D<double>& points, Array1D<double>& funcvalues, Array1D<double>& dervalues);
}

#include <interval/jl_utils.cpp>

#endif

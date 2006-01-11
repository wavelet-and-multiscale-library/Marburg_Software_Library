// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_I_INDEXPLOT_H
#define _WAVELETTL_I_INDEXPLOT_H

#include <iostream>
#include <interval/i_index.h>
#include <utils/plot_tools.h>
#include <algebra/infinite_vector.h>

using MathTL::InfiniteVector;
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*!
    Matlab plot of a given finite wavelet coefficient array
    w.r.t. a wavelet basis on the interval.
    The coefficients should be in multiscale representation, i.e.,
    generator indices are only active on the coarsest level j0.
  */
  template <class IBASIS>
  void plot_indices(const IBASIS* basis,
		    const InfiniteVector<double, typename IBASIS::Index>& coeffs,
		    const int jmax,
		    std::ostream& os,
		    MathTL::MatlabColorMap colormap = jet);
}

// include implementation
#include <interval/i_indexplot.cpp>

#endif

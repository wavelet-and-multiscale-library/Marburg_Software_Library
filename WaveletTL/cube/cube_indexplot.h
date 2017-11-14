// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_INDEXPLOT_H
#define _WAVELETTL_CUBE_INDEXPLOT_H

#include <iostream>
#include <cube/cube_index.h>
#include <algebra/infinite_vector.h>

using MathTL::InfiniteVector;
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*!
    Matlab plot of a given finite wavelet coefficient array
    w.r.t. a wavelet basis on the square.
    The coefficients should be in multiscale representation, i.e.,
    generator indices are only active on the coarsest level j0.

    Since the color value of the rectangles will actually be
    log10(|c_lambda|/||c||_infty), you have to specify a range 10^a...10^0
    to which these values shall be clipped.
    You can choose a Matlab colormap and toggle the boxes and a colorbar on or off.
  */
  template <class CUBEBASIS>
  void plot_indices_cube(const CUBEBASIS* basis,
		    const InfiniteVector<double, typename CUBEBASIS::Index>& coeffs,
		    const int jmax,
		    std::ostream& os,
		    const char* colormap = "cool",
		    bool boxed = false,
		    bool colorbar = true,
		    const double aa = -6);
}

// include implementation
#include <cube/cube_indexplot.cpp>

#endif

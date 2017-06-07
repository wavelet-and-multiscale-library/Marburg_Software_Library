// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TFRAME_INDEXPLOT_H
#define	_WAVELETTL_TFRAME_INDEXPLOT_H

#include <iostream>
#include <cube/tframe_index.h>
#include <algebra/infinite_vector.h>

using MathTL::InfiniteVector;
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*
    Matlab plot of a given finite wavelet coefficient array
    w.r.t. a wavelet frame on the square as modeled by tframe.
    Only indices with 1-norm less or equal than ||j0()||_1+range will be plotted.
    Method works for DIM=2, but the code is already somewhat generalized for
    other dimensions.

    Since the color value of the rectangles will actually be
    log10(|c_lambda|/||c||_infty), you have to specify a range 10^a...10^0
    to which these values shall be clipped.
    You can choose a Matlab colormap and toggle the boxes and a colorbar on or off.
  */
  template <class TENSORFRAME>
  void plot_indices(const TENSORFRAME* frame,
		    const InfiniteVector<double, typename TENSORFRAME::Index>& coeffs,
		    const int range,
		    std::ostream& os,
                    const typename TENSORFRAME::Index::polynomial_type p = typename TENSORFRAME::Index::polynomial_type(),
		    const char* colormap = "cool",
		    bool boxed = false,
		    bool colorbar = true,
		    const double aa = -6);
}

// include implementation
#include <cube/tframe_indexplot.cpp>

#endif	/* _WAVELETTL_TFRAME_INDEXPLOT_H */


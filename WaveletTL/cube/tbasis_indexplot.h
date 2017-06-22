// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TBASIS_INDEXPLOT_H
#define	_WAVELETTL_TBASIS_INDEXPLOT_H

#include <iostream>
#include <cube/tbasis_index.h>
#include <algebra/infinite_vector.h>

using MathTL::InfiniteVector;
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*
    Matlab plot of a given finite wavelet coefficient array
    w.r.t. a wavelet basis on the square as modeled by tbasis.
    Only indices with 1-norm less or equal than ||j0()||_1+range will be plotted.
    Method works for DIM=2, but the code is already somewhat generalized for
    other dimensions.

    Since the color value of the rectangles will actually be
    log10(|c_lambda|/||c||_infty), you have to specify a range 10^a...10^0
    to which these values shall be clipped.
    You can choose a Matlab colormap and toggle the boxes and a colorbar on or off.
  */
  template <class TENSORBASIS>
  void plot_indices(const TENSORBASIS* basis,
		    const InfiniteVector<double, typename TENSORBASIS::Index>& coeffs,
		    const int range,
		    std::ostream& os,
		    const char* colormap = "cool",
		    bool boxed = false,
		    bool colorbar = true,
		    const double aa = -6);
  
  template <class TENSORBASIS>
  void plot_indices2(const TENSORBASIS* basis,
		    const InfiniteVector<double, typename TENSORBASIS::Index>& coeffs,
		    std::ostream& os,
                    const typename TENSORBASIS::Index::level_type& j,
                    const typename TENSORBASIS::Index::type_type& e,
		    const char* colormap = "flipud(gray)",
		    bool boxed = false,
		    bool colorbar = true,
		    const double lowerclim = -6);
}

// include implementation
#include <cube/tbasis_indexplot.cpp>

#endif	/* _WAVELETTL_TBASIS_INDEXPLOT_H */


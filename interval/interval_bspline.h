// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_INTERVAL_BSPLINE_H
#define _WAVELETTL_INTERVAL_BSPLINE_H

#include <numerics/schoenberg_splines.h>
#include <Rd/cdf_utils.h>

using namespace MathTL;

namespace WaveletTL
{
  /*!
    A translated and dilated Schoenberg B-spline as a function object,
    depending on the current value of k, the function is reflected
    at x=1/2, just as in the primal generators of PBasis.
  */
  template <int d>
  class SchoenbergIntervalBSpline_td : public Function<1>
  {
  public:
    //! constructor from j, k
    SchoenbergIntervalBSpline_td(const int j, const int k)
      : j_(j), k_(k) {}
    
    //! virtual destructor
    virtual ~SchoenbergIntervalBSpline_td() {}
    
    //! point value
    inline double value(const Point<1>& p,
			const unsigned int component = 0) const
    {
      if (k_ > (1<<j_)-ell1<d>()-d)
	return EvaluateSchoenbergBSpline_td<d>(j_, (1<<j_)-d-k_-2*ell1<d>(), 1-p[0]);
      else
	return EvaluateSchoenbergBSpline_td<d>(j_, k_, p[0]);
    }
    
    //! point value
    void vector_value(const Point<1> &p,
		      Vector<double>& values) const
    {
      values.resize(1, false);
      values[0] = value(p);
    }
    
    //! set j
    void set_j(const int j) { j_ = j; }
    
    //! set k
    void set_k(const int k) { k_ = k; }
    
  protected:
    int j_, k_;
  };
}

#endif

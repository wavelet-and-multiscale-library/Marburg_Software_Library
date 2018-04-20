// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2018                                            |
// | Philipp Keding, Alexander Sieber                                   |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_SLITDOMAIN_FRAME_EVALUATE_H
#define _WAVELETTL_SLITDOMAIN_FRAME_EVALUATE_H

#include <utils/array1d.h>
#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <slitdomain/slitdomain_frame.h>
//#include <interval/spline_frame.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <class IFRAME> class SlitDomainFrame;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of the slit domain

    We redirect this to a member function of frame.
  */
  template <class IFRAME>
  Array1D<SampledMapping<2> >
  evaluate(const SlitDomainFrame<IFRAME>& frame,
	   const typename SlitDomainFrame<IFRAME>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    return frame.evaluate(lambda, primal, resolution);
  }
   
  /*!
    Evaluate an arbitrary linear combination of primal/dual wavelets
    on a dyadic subgrid of the slit domain

    We redirect this to a member function of frame.
  */
  template <class IFRAME>
  Array1D<SampledMapping<2> >
  evaluate(const SlitDomainFrame<IFRAME>& frame,
	   const InfiniteVector<double, typename SlitDomainFrame<IFRAME>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    return frame.evaluate(coeffs, primal, resolution);
  }
}

#endif


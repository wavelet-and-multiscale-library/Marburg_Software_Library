// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DOMAIN_FRAME_EVALUATE_H
#define _WAVELETTL_DOMAIN_FRAME_EVALUATE_H

#include <utils/array1d.h>
#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <general_domain/domain_frame.h>
//#include <interval/spline_frame.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;

namespace WaveletTL
{
  template <class IFRAME, int NPATCHES> class DomainFrame;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of the L-shaped domain

    We redirect this to a member function of frame.
  */
  template <class IFRAME, int NPATCHES>
  Array1D<SampledMapping<2> >
  evaluate(const DomainFrame<IFRAME, NPATCHES>& frame,
	   const typename DomainFrame<IFRAME, NPATCHES>::Index& lambda,
	   const bool primal,
	   const int resolution)
  {
    return frame.evaluate(lambda, resolution);
  }
   
  /*!
    Evaluate an arbitrary linear combination of primal/dual wavelets
    on a dyadic subgrid of the L-shaped domain

    We redirect this to a member function of frame.
  */
  template <class IFRAME, int NPATCHES>
  Array1D<SampledMapping<2> >
  evaluate(const DomainFrame<IFRAME, NPATCHES>& frame,
	   const InfiniteVector<double, typename DomainFrame<IFRAME, NPATCHES>::Index>& coeffs,
	   const bool primal,
	   const int resolution)
  {
    return frame.evaluate(coeffs, resolution);
  }
}

#endif

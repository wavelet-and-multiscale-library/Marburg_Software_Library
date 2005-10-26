// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_FRAME_EVALUATE_H
#define _FRAMETL_FRAME_EVALUATE_H

#include <geometry/sampled_mapping.h>
#include <algebra/infinite_vector.h>
#include <aggregated_frame.h>

using MathTL::SampledMapping;
using MathTL::InfiniteVector;
using FrameTL::AggregatedFrame;

namespace FrameTL
{
  //template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m> class AggregatedFrame;

  /*!
    Evaluate a single primal/dual generator or wavelet \psi_\lambda
    on a dyadic subgrid of its corresponding patch.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  SampledMapping<DIM_d> evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
			       const typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index& lambda,
			       const bool primal,
			       const int resolution);

  /*!
    Evaluate an arbitrary linear combination of primal/dual wavelets
    on a dyadic subgrid of each patch.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  Array1D<SampledMapping<DIM_d> > evaluate(const AggregatedFrame<IBASIS,DIM_d,DIM_m>& frame,
					   const InfiniteVector<double,
					   typename AggregatedFrame<IBASIS,DIM_d,DIM_m>::Index>& coeffs,
					   const bool primal,
					   const int resolution);
}

#include <frame_evaluate.cpp>

#endif

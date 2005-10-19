// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_AGGREGATED_FRAME_H
#define _FRAMETL_AGGREGATED_FRAME_H

#include <geometry/chart.h>
#include <cube/mapped_cube_basis.h>
#include <geometry/atlas.h>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>

using WaveletTL::MappedCubeBasis;
using MathTL::Atlas;
using MathTL::Array1D;
using MathTL::FixedArray1D;

using namespace WaveletTL;

namespace FrameTL
{
  /*!
    Class for frames on d-dimensional manifolds in R^m.
    The construction principle of these is the following:
    Write the domain or manifold as an overlapping union
    of subdomains (patch), each of them being the smooth
    parametric image of a reference domain, i.e., the
    d-dimensional hypercube. By lifting a wavelet basis
    on the reference domain to the subdomains and taking
    the union of these lifted bases, a frame is obtained.
    
    The manifold is given by an appropriate Atlas. The
    corresponding reference bases, respectively their
    lifted versions, are then be internally constructed.
    For each lifted cube, the user may specify
    2*d boundary conditions.

  */
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class AggregatedFrame
  {
  public:    

    /*!
      destructor
     */
    ~AggregatedFrame();

    /*!
      constructor
     */
    AggregatedFrame(const Atlas<DIM_d, DIM_m>*,
		    const Array1D<FixedArray1D<int,2*DIM_d> >&,
		    const Array1D<FixedArray1D<int,2*DIM_d> >&);

    /*!
      access to the local bases on the i-ht patch
     */
    const MappedCubeBasis<IBASIS, DIM_d, DIM_m>*
    get_local_basis(const unsigned int i);

  protected:
    //! pointer to the underlying atlas
    const Atlas<DIM_d, DIM_m>* atlas_;

    //! boundary conditions
    Array1D<FixedArray1D<int,2*DIM_d> > bc_;

  private:
    Array1D<MappedCubeBasis<IBASIS, DIM_d, DIM_m>* > lifted_bases;

  };

  /*!
    stream output of an AggregatedFrame
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  std::ostream& operator << (std::ostream&,
			     const AggregatedFrame<IBASIS, DIM_d, DIM_m>&);

}

#include "aggregated_frame.cpp"

#endif

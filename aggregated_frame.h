// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_AGGREGATED_FRAME_H
#define _FRAMETL_AGGREGATED_FRAME_H

#include <list>

#include <geometry/chart.h>
#include <cube/mapped_cube_basis.h>
#include <geometry/atlas.h>
#include <utils/fixed_array1d.h>
#include <utils/array1d.h>
#include <frame_index.h>

using std::list;

using WaveletTL::MappedCubeBasis;
using MathTL::Atlas;
using MathTL::Array1D;
using MathTL::FixedArray1D;

//using namespace WaveletTL;

namespace FrameTL
{
  //forward declaration of class AggregatedFrame
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class FrameIndex;

  /*!
    Class for frames on d-dimensional manifolds in R^m.
    The construction principle of these is the following:
    Write the domain or manifold as an overlapping union
    of subdomains (patches), each of them being the smooth
    parametric image of a reference domain, i.e., the
    d-dimensional hypercube. By lifting a wavelet basis
    on the reference domain to the subdomains and taking
    the union of these lifted bases, a frame is obtained.
    
    The manifold is given by an appropriate Atlas. The
    corresponding reference bases, respectively their
    lifted versions, are then internally constructed.
    For each lifted cube, the user may specify
    2*d boundary conditions.

  */
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m = DIM_d>
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
      frame index class
    */
    typedef FrameIndex<IBASIS, DIM_d, DIM_m> Index;

    /*!
      access to the local bases
     */
    const Array1D<MappedCubeBasis<IBASIS, DIM_d, DIM_m>* > bases() const
    {
      return lifted_bases;
    }

    /*!
      geometric type of the support sets.
      support is only implicitly given by the
      support of the corresponding wavelet on the
      hypercube and the respective patch chart.
    */
    typedef struct {
      int j;
      int a[DIM];
      int b[DIM];
      const Chart<DIM_d,DIM_m>* ch;
    } Support;

    /*!
      coarsest level
     */
    const int j0() const { return j0_; }

    /*!
      number of patches
     */
    const unsigned int n_p() const { return lifted_bases.size(); }

    /*!
      access to underlying atlas
     */
    const Atlas<DIM_d, DIM_m>* atlas() const { return atlas_;  }

    /*!
      critical Sobolev regularity for the primal generators/wavelets
    */
    static double primal_regularity() { return IBASIS::primal_regularity(); }
    
    /*!
      number of vanishing moments for the primal wavelets
    */
    static unsigned int primal_vanishing_moments() { return IBASIS::primal_vanishing_moments(); }

    /*!
      number of vanishing moments for the dual wavelets
    */
    static unsigned int dual_vanishing_moments() { return IBASIS::dual_vanishing_moments(); }

    //! index of first generator on level j >= j0
    Index first_generator(const int j) const;
      
    //! index of last generator on level j >= j0
    Index last_generator(const int j) const;
      
    //! index of first wavelet on level j >= j0
    Index first_wavelet(const int j) const;
      
    //! index of last wavelet on level j >= j0
    Index last_wavelet(const int j) const;

  protected:
    //! pointer to the underlying atlas
    const Atlas<DIM_d, DIM_m>* atlas_;

    //! primal boundary conditions
    Array1D<FixedArray1D<int,2*DIM_d> > bc_;

    //! dual boundary conditions
    Array1D<FixedArray1D<int,2*DIM_d> > bcT_;

    //! coarsest possible level j0
    int j0_;

  private:

    /*!
      collection of mapped cube bases together forming
      the aggregated frame
     */
    Array1D<MappedCubeBasis<IBASIS, DIM_d, DIM_m>* > lifted_bases;

    /*!
      the instances of the mapped cube bases 
    */
    list<MappedCubeBasis<IBASIS, DIM_d, DIM_m>*> bases_infact;


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

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_FRAME_INDEX_H
#define _FRAMETL_FRAME_INDEX_H

#include <iostream>
#include <cube/cube_index.h>
#include "aggregated_frame.h"

using std::cout;
using std::endl;

using WaveletTL::CubeIndex;
using WaveletTL::CubeBasis;


namespace FrameTL
{
  /*!
    This class represents (tensor product) multilevel indices.
    The spatial dimension is passed as template parameter. Such a frame index
    mainly consists of a CubeIndex and the respective patchnumber.
   */
  template <class IBASIS, unsigned int DIM_d, unsigend int DIM_m>
  class FrameIndex
  {

  public:
    //! type index type
    typedef MultiIndex<unsigned int,DIM> type_type;
    
    //! translation index type
    typedef MultiIndex<int,DIM> translation_type;

    /*!
      default constructor
    */
    FrameIndex();

    /*!
      copy constructor
    */
    FrameIndex(const FrameIndex&);

    /*!
      constructor
    */
    FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
	       const CubeIndex<IBASIS,DIM>&, const int patch,
	       const unsigned int num_patches_);

    /*!
      constructor
    */
    FrameIndex(const int j,
	       const type_type& e,
	       const translation_type& k,
	       const unsigned int patch,
	       CubeBasis<IBASIS,DIM_d>* basis,
	       const unsigned int num_patches_,
	       const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame);

    
    //! check equality
    bool operator == (const FrameIndex& lambda) const;
    
    //! check non-equality
    inline bool operator != (const FrameIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    FrameIndex& operator ++ ();

    //! lexicographic order <
    bool operator < (const FrameIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const FrameIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    const int j() const { return cbI_.j(); }

    //! type e
    const type_type& e() const { return cbI_.e(); }

    //! translation index k
    const translation_type& k() const { return cbI_.k(); }

    //! access to patchnumber
    const unsigned int p() const { return p_; }

    //! access to CubeIndex
    const CubeIndex<IBASIS,DIM> get_CubeIndex() const
    { return cbI_; }

    //! access to underlying frame
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* get_frame() const
    { return frame_; }

  protected:
    /*!
      patchnumber
     */
    unsigned int p_;

    /*!
      corresponding cube index
     */
    CubeIndex<IBASIS, DIM_d> cbI_;

    /*!
      pointer to corresponding frame
     */
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame_;

  private:
    //total number of patches
    unsigned int num_patches;
  };

  /*!
    stream output for FrameIndex
   */
  template <class IBASIS, unsigned int DIM>
  std::ostream& operator << (std::ostream&, const FrameIndex<IBASIS,DIM>&);


}

// include implementation
#include "frame_index.cpp"

#endif

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
  //forward declaration of class AggregatedFrame
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class AggregatedFrame;

  /*!
    This class represents (tensor product) multilevel indices.
    The spatial dimension is passed as template parameter. Such a frame index
    mainly consists of a CubeIndex and the respective patchnumber.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class FrameIndex
  {

  public:
    //! type index type
    typedef MultiIndex<unsigned int,DIM_d> type_type;
    
    //! translation index type
    typedef MultiIndex<int,DIM_d> translation_type;

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
	       const CubeIndex<IBASIS,DIM_d>&, const int patch);

    /*!
      constructor
    */
    FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
	       const int j,
	       const type_type& e,
	       const unsigned int patch,	       
	       const translation_type& k);
    
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
    const int j() const { return j_; }

    //! type e
    const type_type& e() const { return e_; }

    //! translation index k
    const translation_type& k() const { return k_; }

    //! access to patchnumber
    const unsigned int p() const { return p_; }

    //! access to underlying frame
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* get_frame() const
    { return frame_; }

  protected:

    /*!
      pointer to corresponding frame
     */
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame_;

    //! scale
    int j_;
    
    //! type
    MultiIndex<unsigned int,DIM_d> e_;

    /*!
      patchnumber
     */
    unsigned int p_;

    //! translation
    MultiIndex<int,DIM_d> k_;
    
  };

  /*!
    stream output for FrameIndex
   */
  template <class IBASIS, unsigned int DIM>
  std::ostream& operator << (std::ostream&, const FrameIndex<IBASIS, DIM_d, DIM_m>&);


}

// include implementation
#include "frame_index.cpp"

#endif

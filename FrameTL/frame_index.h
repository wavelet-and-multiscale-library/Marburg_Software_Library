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
#include <aggregated_frame.h>

using std::cout;
using std::endl;

using WaveletTL::CubeIndex;

#include <utils/multiindex.h>

using MathTL::MultiIndex;

namespace FrameTL
{
  
  /*!
    A simple type representing pairs of a frame index, given by
    its number, and a double valued coefficient in a frame expansion.
    It shall be used as replacement for a pair of frame_index and
    coeffient value.
   */
  typedef struct {
    
    //! Will be used to represent the number of a frame index.
    int num;

    //! One coefficient of a frame expansion corresponding to the frame index given by 'num'.
    double val;
    
  } Coefficient;
  

  // forward declaration of class AggregatedFrame
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class AggregatedFrame;

  /*!
    This class represents tensor product multilevel aggregated wavelet frame
    indices. A FrameIndex mainly encodes level, type, patchnumber and shift parameter of an
    aggregated wavelet frame element.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m = DIM_d>
  class FrameIndex
  {

  public:
    
    //! The index type.
    typedef MultiIndex<int,DIM_d> type_type;
    
    //! The translation index type.
    typedef MultiIndex<int,DIM_d> translation_type;

    /*!
      Default constructor.
    */
    FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame = 0);

    /*!
      Copy constructor.
    */
    FrameIndex(const FrameIndex&);

    /*!
      Constructor. Clone index from a given const pointer to another FrameIndex.
     */
    FrameIndex(const FrameIndex*);
    
    /*!
      Constructor. Creates the FrameIndex from a given wavelet index of
      a wavelet on the unit cube. Also the integer number (according to the
      canonical lexicographical ordering "level, type, patch, shift parameter")
      of this FrameIndex is calculated.
    */
    FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
	       const CubeIndex<IBASIS,DIM_d>&, const int patch);

    /*!
      Constructor. Create a FrameIndex from the given level, type patchnumber
      and shift parameter. Also the integer number (according to the
      canonical lexicographical ordering "level, type, patch, shift parameter")
      of this FrameIndex is calculated.
    */
    FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
	       const int j,
	       const type_type& e,
	       const unsigned int patch,	       
	       const translation_type& k);

    /*!
      Constructor. From the given integer number (according to the
      canonical lexicographical ordering "level, type, patch, shift parameter")
      the a FrameIndex is reconstructed, i.e., level, type, patchnumber, and
      shift parameter are determined.
    */
    FrameIndex(const int num,
	       const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame);

    
    //! Check for equality.
    bool operator == (const FrameIndex& lambda) const;
    
    //! Check for non-equality.
    inline bool operator != (const FrameIndex& lambda) const
    { return !(*this == lambda); }
    
    //! Preincrement.
    FrameIndex& operator ++ ();

    //! Assignment.
    FrameIndex& operator = (const FrameIndex&);

    //! Lexicographic order < .
    bool operator < (const FrameIndex& lambda) const;

    //! Lexicographic order <= .
    bool operator <= (const FrameIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! Scale j.
    const int j() const { return j_; }

    //! Type e.
    const type_type& e() const { return e_; }

    //! Translation index k.
    const translation_type& k() const { return k_; }

    //! Access to patchnumber.
    const int p() const { return p_; }

    const int number() const { return num_; }

    //! Access to underlying frame.
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame() const
    { return frame_; }

    /*!
      The integer number (according to the
      canonical lexicographical ordering "level, type, patch, shift parameter")
      of this FrameIndex is calculated and set.
    */
    void set_number();

  protected:
    
    /*!
      Pointer to corresponding frame.
     */
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame_;

    //! Scale.
    int j_;
    
    //! Type.
    type_type e_;

    /*!
      Patchnumber.
     */
    int p_;

    //! Translation index.
    translation_type k_;

    //! Number of this index.
    int num_;

  };

  /*!
    Stream output for FrameIndex.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  std::ostream& operator << (std::ostream&, const FrameIndex<IBASIS, DIM_d, DIM_m>&);

  /*!
    Index of first generator on level \f$j \geq j_0\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  first_generator(const FRAMETYPE* frame, const int j);

  /*!
    Index of last generator on level \f$j \geq j_0\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  last_generator(const FRAMETYPE* frame,const int j);

    
  /*!
    Index of first wavelet on level \f$j \geq j_0\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  first_wavelet(const FRAMETYPE* frame, const int j);
    
  /*!
    Index of last wavelet on level \f$j \geq j_0\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  last_wavelet(const FRAMETYPE* frame, const int j);

  /*!
    Number of first generator on level \f$j_0\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  int
  first_generator_num(const FRAMETYPE* frame);

  /*!
    Number of last generator on level \f$j_0\f$.
    */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  int
  last_generator_num(const FRAMETYPE* frame);
    
  /*!
    Number of first wavelet on level \f$j \geq j0\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  int
  first_wavelet_num(const FRAMETYPE* frame, const int j);
    
  /*!
    Index of last wavelet on level \f$j \geq j0\f$.
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAMETYPE>
  int
  last_wavelet_num(const FRAMETYPE* frame, const int j);
  

  // The following are helpers for the parallel implementation, where,
  // to send an InfiniteVector from one processor to another one,
  // the InfiniteVector is first copied into an array of Coefficient's
  // which is then sent to the recipient and copied there into an
  // InfiniteVector again.
  template <class INDEX>
  void to_array (const InfiniteVector<double, INDEX>& ivec,
		 Coefficient* coeff_array);
  
  template <class INDEX, class FRAMETYPE>
  void array_to_map (const Coefficient* coeff_array,
                     const FRAMETYPE* frame,
		     InfiniteVector<double, INDEX>& ivec,
		     const int count);

}
// include implementation
#include "frame_index.cpp"

#endif

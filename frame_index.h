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
    It shall be used as repacement for a pair of frame_index and
    coeffient value.
   */
  typedef struct {
    
    // will be used to represent the number of a frame index
    int num;

    // one coefficient of a frame expansion corresponding to the
    // frame index given by 'num'
    double val;
    
  } Coefficient;
  

  //forward declaration of class AggregatedFrame
  template<class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  class AggregatedFrame;

  /*!
    This class represents (tensor product) multilevel indices.
    The spatial dimension is passed as template parameter. Such a frame index
    mainly consists of a CubeIndex and the respective patchnumber.
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m = DIM_d>
  class FrameIndex
  {

  public:
    //! type index type
    typedef MultiIndex<int,DIM_d> type_type;
    
    //! translation index type
    typedef MultiIndex<int,DIM_d> translation_type;

    /*!
      default constructor
    */
    FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame = 0);

    /*!
      copy constructor
    */
    FrameIndex(const FrameIndex&);

    /*!
      clone index from a given const pointer to another FrameIndex
     */
    FrameIndex(const FrameIndex*);


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

    /*!
      constructor
    */
    FrameIndex(const int num,
	       const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame);

    
    //! check equality
    bool operator == (const FrameIndex& lambda) const;
    
    //! check non-equality
    inline bool operator != (const FrameIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    FrameIndex& operator ++ ();

    //! assignment
    FrameIndex& operator = (const FrameIndex&);

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
    const int p() const { return p_; }

    const int number() const { return num_; }

    //! access to underlying frame
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame() const
    { return frame_; }

    /*!
      inverse of constructor
      'FrameIndex(const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame,
      const unsigned int num)'
    */
    void set_number();

  protected:
    
    /*!
      pointer to corresponding frame
     */
    const AggregatedFrame<IBASIS,DIM_d,DIM_m>* frame_;

    //! scale
    int j_;
    
    //! type
    //MultiIndex<unsigned int,DIM_d> e_;
    type_type e_;

    /*!
      patchnumber
     */
    int p_;

    //! translation
    //MultiIndex<int,DIM_d> k_;
    translation_type k_;
    
    int num_;

  };

  /*!
    stream output for FrameIndex
   */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m>
  std::ostream& operator << (std::ostream&, const FrameIndex<IBASIS, DIM_d, DIM_m>&);

  /*!
    index of first generator on level j >= j0
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  first_generator(const FRAME* frame, const int j);

  /*!
    index of last generator on level j >= j0
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  last_generator(const FRAME* frame,const int j);

    
  /*!
    index of first wavelet on level j >= j0
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  first_wavelet(const FRAME* frame, const int j);
    
  /*!
    index of last wavelet on level j >= j0
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  FrameIndex<IBASIS,DIM_d,DIM_m>
  last_wavelet(const FRAME* frame, const int j);

  /*!
    number of first generator on level j0
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  int
  first_generator_num(const FRAME* frame);

  /*!
    number of last generator on level j0
    */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  int
  last_generator_num(const FRAME* frame);
    
  /*!
    number of first wavelet on level j >= j0
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  int
  first_wavelet_num(const FRAME* frame, const int j);
    
  /*!
    index of last wavelet on level j >= j0
  */
  template <class IBASIS, unsigned int DIM_d, unsigned int DIM_m, class FRAME>
  int
  last_wavelet_num(const FRAME* frame, const int j);
  
  template <class INDEX>
  void to_array (const InfiniteVector<double, INDEX>& ivec,
		 Coefficient* coeff_array);
  
  template <class INDEX, class FRAME>
  void array_to_map (const Coefficient* coeff_array,
		     const FRAME* frame,
		     InfiniteVector<double, INDEX>& ivec,
		     const int count);

}
// include implementation
#include "frame_index.cpp"

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_HAAR_MASK_H
#define _WAVELETTL_HAAR_MASK_H

#include <Rd/r_mask.h>

namespace WaveletTL
{
  class HaarMask
    : public virtual RRefinementMask<2>
  {
  public:
    //! default constructor, sets up the mask
    HaarMask() { RRefinementMask<2>::coeffs[0] = RRefinementMask<2>::coeffs[1] = 1.0; }
    
    //! Strang-Fix order, i.e., degree of polynomial reproduction (+1)
    static unsigned int Strang_Fix_order() { return 1; }
    
    //! critical L_2 Sobolev regularity
    static double regularity() { return 0.5; }  

    //! index of first nontrivial refinement coefficient (inherited from RRefinementMask)
    const int abegin() const { return 0; }
  };
}

#endif

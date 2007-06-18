// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_INDEX1D_H
#define _FRAMETL_INDEX1D_H

namespace FrameTL
{
  template <class IBASIS>
  class Index1D
  {

  public:
    Index1D (const typename IBASIS::Index& ind,
	     const unsigned int p, const unsigned int dir,
	     const unsigned int der);

    bool operator < (const Index1D<IBASIS>& lambda) const;
    bool operator == (const Index1D<IBASIS>& lambda) const;
    bool operator != (const Index1D<IBASIS>& lambda) const;
    bool operator <= (const Index1D<IBASIS>& lambda) const;

    typename IBASIS::Index index() const { return ind_; };
    unsigned int derivative() const { return der_; };
    unsigned int p() const { return p_; };
    unsigned int direction() const { return dir_; };

  protected:
    typename IBASIS::Index ind_;
    unsigned int p_;
    unsigned int dir_;
    unsigned int der_;

  };
}
#include <index1D.cpp>

#endif

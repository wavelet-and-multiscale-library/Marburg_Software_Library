// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_R_MASK_H
#define _WAVELETTL_R_MASK_H

#include <iostream>
#include <utils/fixed_array1d.h>
#include <geometry/sampled_mapping.h>

using MathTL::FixedArray1D;
using MathTL::SampledMapping;

namespace WaveletTL
{
  /*!
    base class for 1D refinement masks of length L
  */
  template <unsigned int L>
  class RRefinementMask;

  //! stream output for a refinement mask
  template <unsigned int L>
  std::ostream& operator << (std::ostream& os, const RRefinementMask<L>& m)
  {
    os << m.coeffs;
    return os;
  }

  template <unsigned int L>
  class RRefinementMask
  {
    friend std::ostream& operator << <L>(std::ostream& os, const RRefinementMask<L>& m);

  public:
    //! make template argument accessible
    static const unsigned int length = L;
    
    //! virtual destructor
    virtual ~RRefinementMask() {}
    
    //! index of first nontrivial refinement coefficient
    virtual const int abegin() const = 0;

    //! index of last nontrivial refinement coefficient
    const int aend() const { return abegin()+L-1; }
    
    //! k-th refinement coefficient
    const double a(const int k) const { return (k<abegin() ? 0 : (k>aend() ? 0 : coeffs[k-abegin()])); }

  protected:
    //! refinement coefficients
    FixedArray1D<double,L> coeffs;
  };
  
}

#endif

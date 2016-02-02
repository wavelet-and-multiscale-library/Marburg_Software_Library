// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_R_MASK_H
#define _WAVELETTL_R_MASK_H

#include <iostream>
#include <utils/fixed_array1d.h>
#include <algebra/infinite_vector.h>
#include <geometry/sampled_mapping.h>

using MathTL::FixedArray1D;
using MathTL::InfiniteVector;
using MathTL::SampledMapping;

namespace WaveletTL
{
  /*!
    base class for 1D refinement masks of length L and offset BEGIN
  */
  template <unsigned int L, int BEGIN>
  class RRefinementMask;

  //! stream output for a refinement mask
  template <unsigned int L, int BEGIN>
  std::ostream& operator << (std::ostream& os, const RRefinementMask<L, BEGIN>& m)
  {
    os << m.coeffs_;
    return os;
  }

  template <unsigned int L, int BEGIN>
  class RRefinementMask
  {
    friend std::ostream& operator << <L, BEGIN> (std::ostream& os, const RRefinementMask<L, BEGIN>& m);

  public:
    //! make template argument accessible
    static const unsigned int length = L;
    
    //! virtual destructor
    virtual ~RRefinementMask() {}
    
    //! index of first nontrivial refinement coefficient
    static int begin() { return BEGIN; }

    //! index of last nontrivial refinement coefficient
    static int end() { return begin()+L-1; }
    
    //! k-th refinement coefficient
    double a(const int k) const { return (k<begin() ? 0 : (k>end() ? 0 : coeffs_[k-begin()])); }
    
    //! read access to the full coefficient array
    const FixedArray1D<double,L>& coeffs() const { return coeffs_; }
    
    /*!
      Evaluate the refinable function \phi on the dyadic grid 2^{-resolution}\mathbb Z.
      We assume that \phi is zero at the boundary of its support.
    */
    InfiniteVector<double, int>
    evaluate(const int resolution = 0) const;

    /*!
      Evaluate the mu-th derivative of the refinable function \phi
      on the dyadic grid 2^{-resolution}\mathbb Z.
      We assume that the derivative of \phi is zero at the boundary of its support (!).
    */
    InfiniteVector<double, int>
    evaluate(const int mu, const int resolution) const;
    
    /*!
      Evaluate the mu-th derivative of a dilated and translated version
      of the refinable function \phi
      
      (d/dx)^\mu \phi_{j,k}(x) = 2^{j\mu} * 2^{jd/2} * \phi^{(\mu)}(2^j*x-k)
      
      on a dyadic subgrid of the interval [a,b].
      
      We assume that \partial^\mu\phi is zero at the boundary of its support.
    */
    SampledMapping<1>
    evaluate(const int mu,
	     const int j,
	     const int k,
	     const int a,
	     const int b,
	     const int resolution) const;

    /*!
      Compute the k-th moment of the refinable function \phi via the mask
    */
    const double moment(const unsigned int k) const;

  protected:
    //! refinement coefficients
    FixedArray1D<double,L> coeffs_;
  };
  
}

#include <Rd/r_mask.cpp>

#endif

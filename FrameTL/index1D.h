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

  /*!
    This class is just used to attach to an IntervalIndex<IBASIS> from interval/i_index.h
    a patchnumber, a spatial direction and a derivative order. This class is only 
    used, e.g., in simple_elliptic_equation.cpp (routine integrate) for the caching of
    one dimensional integrals (in a certain spatial direction)
    arising when the tensor product structure is used to
    evaluate the bilinear form associated to an elliptic operator.
   */
  template <class IBASIS>
  class Index1D
  {
  public:
    
    /*!
      Constructor.
     */
    Index1D (const typename IBASIS::Index& ind,
	     const unsigned int p, const unsigned int dir,
	     const unsigned int der);

    /*!
      Relation "less than" according to the hierarchy
      "IntervalIndex, patch, direction, derivative".
     */
    bool operator < (const Index1D<IBASIS>& lambda) const;

    /*!
      Check for equality.
     */
    bool operator == (const Index1D<IBASIS>& lambda) const;

    /*!
      Check for inequality.
     */
    bool operator != (const Index1D<IBASIS>& lambda) const;

    /*!
      Relation "less than or equal" according to the hierarchy
      "IntervalIndex, patch, direction, derivative".
     */
    bool operator <= (const Index1D<IBASIS>& lambda) const;

    /*!
      Read access to the encapsulated 1D wavelet index.
     */
    typename IBASIS::Index index() const { return ind_; };

    /*!
      The order of derivative of the encoded wavelet or generator.
     */
    unsigned int derivative() const { return der_; };
    
    
    /*!
      The patchnumber of the encoded wavelet or generator.
     */
    unsigned int p() const { return p_; };
   
    /*!
      The spatial direction the 1D wavelet or generator belongs to.
     */
    unsigned int direction() const { return dir_; };

  protected:

    /*!
      The 1D wavelet or generator index.
     */
    typename IBASIS::Index ind_;
    
    /*!
      The patchnumber.
     */
    unsigned int p_;

    /*!
      The spatial direction.
     */
    unsigned int dir_;

    /*!
      The derivative order.
     */
    unsigned int der_;

  };
}
#include <index1D.cpp>

#endif

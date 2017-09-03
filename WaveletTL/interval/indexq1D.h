// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner, Philipp Keding                                      |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_INDEXQ1D_H
#define _WAVELETTL_INDEXQ1D_H

namespace WaveletTL
{

  /*!
    This class is just used to attach to an IntervalQIndex<IFRAME> from interval/i_q_index.h
    a patchnumber, a spatial direction and a derivative order. This class is only 
    used, e.g., in ldomain_frame_equation.cpp (routine integrate) for the caching of
    one dimensional integrals (in a certain spatial direction)
    arising when the tensor product structure is used to
    evaluate the bilinear form associated to an elliptic operator.
   */
  template <class IFRAME>
  class IndexQ1D
  {
  public:
    
    /*!
      Constructor.
     */
    IndexQ1D (const typename IFRAME::Index& ind,
	     const unsigned int patch, const unsigned int dir,
	     const unsigned int der);

    /*!
      Relation "less than" according to the hierarchy
      "IntervalIndex, patch, direction, derivative".
     */
    bool operator < (const IndexQ1D<IFRAME>& lambda) const;

    /*!
      Check for equality.
     */
    bool operator == (const IndexQ1D<IFRAME>& lambda) const;

    /*!
      Check for inequality.
     */
    bool operator != (const IndexQ1D<IFRAME>& lambda) const;

    /*!
      Relation "less than or equal" according to the hierarchy
      "IntervalIndex, patch, direction, derivative".
     */
    bool operator <= (const IndexQ1D<IFRAME>& lambda) const;

    /*!
      Read access to the encapsulated 1D wavelet index.
     */
    typename IFRAME::Index index() const { return ind_; };

    /*!
      The order of derivative of the encoded wavelet or generator.
     */
    unsigned int derivative() const { return der_; };
    
    
    /*!
      The patchnumber of the encoded wavelet or generator.
     */
    unsigned int patch() const { return patch_; };
   
    /*!
      The spatial direction the 1D wavelet or generator belongs to.
     */
    unsigned int direction() const { return dir_; };

  protected:

    /*!
      The 1D wavelet or generator index.
     */
    typename IFRAME::Index ind_;
    
    /*!
      The patchnumber.
     */
    unsigned int patch_;

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
#include <indexq1D.cpp>

#endif

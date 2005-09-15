// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Manuel Werner                                                      |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_ATLAS_H
#define _FRAMETL_ATLAS_H

#include <iostream>

using std::cout;
using std::endl;

namespace FrameTL
{
  /*!
    This class provides a realization of an arbitrary bounded domain
    in \mathbb R^{DIM_d}, that is covered with (overlapping) patches each of
    which are represented by a mapping from \mathbb R^{DIM_m} onto
    \mathbb R^{DIM_d}.
  */
  template <unsigned int DIM_m, unsigned int DIM_d>
  class Atlas
  {

  public:
    /*!
      default constructor
     */
    Atlas () { };

    /*!
      copy constructor
     */
    Atlas (const Atlas &);

    /*!
      constructor
     */
    Atlas (const Vector<Parametrization<DIM_m, DIM_d> *>& charts);

  protected:
    /*!
      number of patches
     */
    unsigned int num_patches;

    /*!
      the parametrizations of the patches
     */
    Vector<Parametrization<DIM_m, DIM_d> *> charts;

  private:

  };

  /*!
    stream output for Atlas
   */
  template <unsigned int DIM_m, unsigned int DIM_d>
  std::ostream& operator << (std::ostream&, const Atlas<DIM_m, DIM_d>&);

}

// include implementation
#include "atlas.cpp"

#endif


// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_UTILS_H
#define _WAVELETTL_CDF_UTILS_H

namespace WaveletTL
{
  //! left support bound for a primal CDF generator
  template <int d>
    int ell1() { return -(d/2); }
  
  //! right support bound for a primal CDF generator
  template <int d>
    int ell2() { return d - (d/2); }
  
  //! left support bound for a dual CDF generator
  template <int d, int dt>
    int ell1T() { return ell1<d>() - dt + 1; }
  
  //! right support bound for a dual CDF function
  template <int d, int dt>
    int ell2T() { return ell2<d>() + dt - 1; }
}

#endif

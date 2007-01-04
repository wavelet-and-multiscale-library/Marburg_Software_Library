// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CDF_UTILS_H
#define _WAVELETTL_CDF_UTILS_H

#include <Rd/r_index.h>

namespace WaveletTL
{
  //! left support bound for a primal CDF generator
  template <int d>
  int ell1() { return -(d/2); }
  
  //! right support bound for a primal CDF generator
  template <int d>
  int ell2() { return d - (d/2); }
  
  //! left support bound for a dual CDF generator
  template <int d, int dT>
  int ell1T() { return ell1<d>() - dT + 1; }
  
  //! right support bound for a dual CDF generator
  template <int d, int dT>
  int ell2T() { return ell2<d>() + dT - 1; }

  //! left support bound for a primal (or dual) CDF wavelet
  template <int d, int dT>
  int psi_supp_left() { return -(d + dT) / 2 + 1; }

  //! right support bound for a primal (or dual) CDF wavelet
  template <int d, int dT>
  int psi_supp_right() { return (d + dT) / 2; }

  /*!
    compute support bounds of the form 2^{-j}k for a (translated and dilated)
    CDF generator or wavelet
    (j == lambda.j() is neglected for performance reasons)
  */
  template <int d, int dT>
  void support(const RIndex& lambda, const bool primal, int& k1, int& k2) {
    k1 = (lambda.e() == 0 ? (primal ? ell1<d>() : ell1T<d,dT>()) : psi_supp_left<d,dT>()) + lambda.k();
    k2 = (lambda.e() == 0 ? (primal ? ell2<d>() : ell2T<d,dT>()) : psi_supp_right<d,dT>()) + lambda.k();
  }
}

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TP_BASIS_H
#define _WAVELETTL_TP_BASIS_H

namespace WaveletTL
{
  /*!
    Template class for a tensor product wavelet basis from two wavelet bases
    Psi, Xi over bounded domains.
    The template parameters BASISi may or may not allow the specification of
    boundary conditions.
  */
  template <class BASIS1, class BASIS2>
  class TPBasis
  {
  public:

  protected:
    // instances of the 1D bases
    BASIS1 basis1;
    BASIS2 basis2;
  };
}

#include <generic/tp_basis.cpp>

#endif

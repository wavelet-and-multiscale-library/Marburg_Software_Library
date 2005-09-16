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
  class TensorProductBasis
  {
  public:
    //! default constructor
    TensorProductBasis();

    //! coarsest possible level j0
    inline const int j0() const { return j0_; }    

  protected:
    //! coarsest possible level j0
    int j0_;

    //! instances of the two 1D bases
    BASIS1 basis1_;
    BASIS2 basis2_;
  };
}

#include <generic/tp_basis.cpp>

#endif

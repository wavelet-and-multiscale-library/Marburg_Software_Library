// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TP_BASIS_H
#define _WAVELETTL_TP_BASIS_H

#include <generic/tp_index.h>

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

    //! wavelet index class
    typedef TensorProductIndex<BASIS1,BASIS2> Index;

    //! read access to the first basis
    const BASIS1& basis1() const { return basis1_; }

    //! read access to the second basis
    const BASIS2& basis2() const { return basis2_; }

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

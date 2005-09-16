// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TP_INDEX_H
#define _WAVELETTL_TP_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*!
    An index class for tensor product wavelet bases over bounded domains.
    The template parameter bases BASISi should provide (via typedef) the
    index classes BASISi::Index, i=1,2.
  */
  template <class BASIS1, class BASIS2>
  class TensorProductIndex
  {
  public:
    /*!
      constructor with given bases
      (also serves as a default constructor, but yields an invalid index then
      because the underlying bases must be specified to work correctly)
     */
    IntervalIndex(const BASIS1* basis1 = 0, const BASIS2* basis2 = 0);

  protected:
    //! pointer to the underlying bases
    const BASIS1* basis1_;
    const BASIS2* basis2_;

    //! current entries
    typename BASIS1::Index index1;
    typename BASIS2::Index index2;
  };
}

// include implementation
#include <generic/tp_index.cpp>

#endif

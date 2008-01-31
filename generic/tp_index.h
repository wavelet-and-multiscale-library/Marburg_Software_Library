// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_TP_INDEX_H
#define _WAVELETTL_TP_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  // forward declaration
  template <class BASIS0, class BASIS1>
  class TensorProductBasis;

  /*!
    An index class for tensor product wavelet bases over bounded domains.
    The template parameter bases BASISi should provide (via typedef) the
    index classes BASISi::Index, i=0,1.
  */
  template <class BASIS0, class BASIS1>
  class TensorProductIndex
  {
  public:
    /*!
      default constructor
     */
    TensorProductIndex();

    //! copy constructor
    TensorProductIndex(const TensorProductIndex& lambda);
  
    //! constructor from a tensor product basis and two single indices
    TensorProductIndex(const typename BASIS0::Index& index0,
		       const typename BASIS1::Index& index1);

    //! assignment
    TensorProductIndex& operator = (const TensorProductIndex& lambda);

    //! check equality
    bool operator == (const TensorProductIndex& lambda) const;
    
    //! check non-equality
    inline bool operator != (const TensorProductIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    TensorProductIndex& operator ++ ();

    //! lexicographic order <
    bool operator < (const TensorProductIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const TensorProductIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    int j() const { return index0_.j(); }

    //! read access to the first index
    const typename BASIS0::Index& index0() const { return index0_; }

    //! read access to the second index
    const typename BASIS1::Index& index1() const { return index1_; }

  protected:
    //! current entries
    typename BASIS0::Index index0_;
    typename BASIS1::Index index1_;
  };

  //! stream output
  template <class BASIS0, class BASIS1>
  inline std::ostream& operator << (std::ostream& os,
				    const TensorProductIndex<BASIS0,BASIS1>& lambda)
  {
    using namespace std;
    os << "("
       << lambda.index0()
       << ","
       << lambda.index1()
       << ")";
    return os;
  }

//   /*!
//     index of first generator on level j >= j0
//   */
//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   first_generator(const int j);

//   /*!
//     index of last generator on level j >= j0
//   */
//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   last_generator(const int j);

//   /*!
//     index of first wavelet on level j >= j0
//   */
//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   first_wavelet(const int j);

//   /*!
//     index of last wavelet on level j >= j0
//   */
//   template <class BASIS0, class BASIS1>
//   TensorProductIndex<BASIS0,BASIS1>
//   last_wavelet(const int j);
}

// include implementation
#include <generic/tp_index.cpp>

#endif

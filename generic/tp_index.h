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
  template <class BASIS1, class BASIS2>
  class TensorProductBasis;

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
      constructor with a given tensor product basis
      (also serves as a default constructor, but yields an invalid index pair
      in this case, because the underlying bases must be specified to work correctly)
     */
    TensorProductIndex(const TensorProductBasis<BASIS1,BASIS2>* basis = 0);

    //! copy constructor
    TensorProductIndex(const TensorProductIndex& lambda);
  
    //! constructor from a tensor product basis and two single indices
    TensorProductIndex(const TensorProductBasis<BASIS1,BASIS2>* basis,
		       const typename BASIS1::Index& index1,
		       const typename BASIS2::Index& index2);

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
    int j() const { return index1_.j(); }

    //! read access to the first index
    const typename BASIS1::Index& index1() const { return index1_; }

    //! read access to the second index
    const typename BASIS2::Index& index2() const { return index2_; }

  protected:
    //! pointer to the underlying basis
    const TensorProductBasis<BASIS1,BASIS2>* basis_;

    //! current entries
    typename BASIS1::Index index1_;
    typename BASIS2::Index index2_;
  };

  //! stream output
  template <class BASIS1, class BASIS2>
  inline std::ostream& operator << (std::ostream& os,
				    const TensorProductIndex<BASIS1,BASIS2>& lambda)
  {
    using namespace std;
    os << "("
       << lambda.index1()
       << ","
       << lambda.index2()
       << ")";
    return os;
  }

  /*!
    index of first generator on level j >= j0
  */
  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  first_generator(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j);

  /*!
    index of last generator on level j >= j0
  */
  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  last_generator(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j);

  /*!
    index of first wavelet on level j >= j0
  */
  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  first_wavelet(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j);

  /*!
    index of last wavelet on level j >= j0
  */
  template <class BASIS1, class BASIS2>
  TensorProductIndex<BASIS1,BASIS2>
  last_wavelet(const TensorProductBasis<BASIS1,BASIS2>* basis, const int j);
}

// include implementation
#include <generic/tp_index.cpp>

#endif

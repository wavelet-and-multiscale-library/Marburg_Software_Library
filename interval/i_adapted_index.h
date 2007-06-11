// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_I_ADAPTED_INDEX_H
#define _WAVELETTL_I_ADAPTED_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*!
    A (template) wrapper class for adapting a multiwavelet index to a
    classical wavelet index.
  */
  template <class IAdaptedBasis>
  class IntervalAdaptedIndex
  {
  typedef typename IAdaptedBasis::MultiIndex IMultiIndex;

  public:
    /*!
      constructor with given adapted interval basis
      (also serves as a default constructor, but yields an invalid index then
      because the underlying interval basis must be specified to work correctly)
    */
    IntervalAdaptedIndex(const IAdaptedBasis* basis = 0);

    //! construct encoded index from underlying multiwavelet index
    IntervalAdaptedIndex(const IMultiIndex& mu, const IAdaptedBasis* basis);
  
    //! copy constructor
    IntervalAdaptedIndex(const IntervalAdaptedIndex& lambda);

    //! copy index from const pointer
    IntervalAdaptedIndex(const IntervalAdaptedIndex* lambda);

    //! constructor with specified parameters
    IntervalAdaptedIndex(const int j, const int e, const int k, const IAdaptedBasis* basis);

    //! constructor with specified number
    IntervalAdaptedIndex(const int num, const IAdaptedBasis* basis);

    //! destructor
    ~IntervalAdaptedIndex();
    
    //! assignment
    IntervalAdaptedIndex& operator = (const IntervalAdaptedIndex& lambda);

    //! check equality
    bool operator == (const IntervalAdaptedIndex& lambda) const;

    //! check non-equality
    bool operator != (const IntervalAdaptedIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    IntervalAdaptedIndex& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const IntervalAdaptedIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const IntervalAdaptedIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    const int j() const { return index_->j(); }

    //! type index type
    typedef int type_type;

    //! type e
    const type_type& e() const { return index_->e(); }

    //! translation index type
    typedef int translation_type;

    //! translation index k
    const translation_type k() const { return index_->ck_encode(); }

    //! number
    const int number() const {
      return index_->number();
    }

    /*!
      inverse of constructor
      'IntervalAdaptedIndex(const unsigned int number, const IAdaptedBasis* basis)'
    */
    void set_number();


    //! get multi index
    const IMultiIndex* multi_index() const { return index_; }

    //! underlying basis
    const IAdaptedBasis* basis() const { return basis_; }

  protected:  

    //! pointer to an object of the adapted multiwavelet basis
    IMultiIndex* index_;

    //! pointer to the underlying interval basis
    const IAdaptedBasis* basis_;
  };

  //! stream output
  template <class IAdaptedBasis>
  std::ostream& operator << (std::ostream& os, const IntervalAdaptedIndex<IAdaptedBasis>& lambda)
  {
    using namespace std;
    if (lambda.multi_index() == 0)
      os << "NULL";
    else
      os << "("
         << lambda.j()
         << ","
         << lambda.e()
         << ","
         << lambda.k()
         << ")";
    return os;
  }

  /*!
    index of first (leftmost) generator on level j >= j0
  */
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis> first_generator(const IAdaptedBasis* basis, const int j);

  /*!
    index of last (rightmost) generator on level j >= j0
  */
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis> last_generator(const IAdaptedBasis* basis, const int j);

  /*!
    index of first (leftmost) wavelet on level j >= j0
  */
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis> first_wavelet(const IAdaptedBasis* basis, const int j);

  /*!
    index of last (rightmost) wavelet on level j >= j0
  */
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis> last_wavelet(const IAdaptedBasis* basis, const int j);

  /*!
    index of first function with type e
    (mainly for TensorProductBasis)
  */
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis> first_index(const IAdaptedBasis* basis, const int j, const int e);

  /*!
    index of last function with type e
    (mainly for TensorProductBasis)
  */
  template <class IAdaptedBasis>
  IntervalAdaptedIndex<IAdaptedBasis> last_index(const IAdaptedBasis* basis, const int j, const int e);

  /*!
    number of first (leftmost) generator on level j0
  */
  template <class IAdaptedBasis>
  int first_generator_num(const IAdaptedBasis* basis);

  /*!
    number of last (rightmost) generator on level j0
  */
  template <class IAdaptedBasis>
  int last_generator_num(const IAdaptedBasis* basis);

  /*!
    number of first (leftmost) wavelet on level j >= j0
  */
  template <class IAdaptedBasis>
  int first_wavelet_num(const IAdaptedBasis* basis, const int j);

  /*!
    number of last (rightmost) wavelet on level j >= j0
  */
  template <class IAdaptedBasis>
  int last_wavelet_num(const IAdaptedBasis* basis, const int j);

}

// include implementation
#include <interval/i_adapted_index.cpp>

#endif // _WAVELETTL_I_ADAPTED_INDEX_H

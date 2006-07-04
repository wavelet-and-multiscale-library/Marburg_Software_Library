// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_I_INDEX_H
#define _WAVELETTL_I_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*!
    A (template) wavelet index class for bases on [0,1] like periodic wavelet bases
    or those from [DKU],[DS] (see template class DSBasis).
    We require the template parameter class IBASIS to provide the routines
    j0(), DeltaLmin(), DeltaRmax(j), Nablamin() and Nablamax(j).
  */
  template <class IBASIS>
  class IntervalIndex
  {
  public:
    /*!
      constructor with given interval basis
      (also serves as a default constructor, but yields an invalid index then
      because the underlying interval basis must be specified to work correctly)
    */
    IntervalIndex(const IBASIS* basis = 0);

    //! copy constructor
    IntervalIndex(const IntervalIndex& lambda);
  
    //! constructor with specified parameters
    IntervalIndex(const int j, const int e, const int k, const IBASIS* basis);

    //! constructor with specified number
    IntervalIndex(const int num, const IBASIS* basis);
    
    //! assignment
    IntervalIndex& operator = (const IntervalIndex& lambda);

    //! check equality
    bool operator == (const IntervalIndex& lambda) const;

    //! check non-equality
    bool operator != (const IntervalIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    IntervalIndex& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const IntervalIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const IntervalIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    const int j() const { return j_; }

    //! type index type
    typedef int type_type;

    //! type e
    const type_type& e() const { return e_; }

    //! translation index type
    typedef int translation_type;

    //! translation index k
    const translation_type& k() const { return k_; }

    //! number num_
    const int number() const { return num_; }

    /*!
      inverse of constructor
      'IntervalIndex(const unsigned int number,
                     const IBASIS* basis)'
    */
    void set_number();


    //! underlying basis
    const IBASIS* basis() const { return basis_; }

  protected:  
    
    //! scale, type, translation
    int j_, e_, k_;
    
    //! number of the index (not for generators on level j > j0)
    int num_;
    

    //! pointer to the underlying interval basis
    const IBASIS* basis_;
  };

  //! stream output
  template <class IBASIS>
  std::ostream& operator << (std::ostream& os, const IntervalIndex<IBASIS>& lambda)
  {
    using namespace std;
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
  template <class IBASIS>
  IntervalIndex<IBASIS> first_generator(const IBASIS* basis, const int j);

  /*!
    index of last (rightmost) generator on level j >= j0
  */
  template <class IBASIS>
  IntervalIndex<IBASIS> last_generator(const IBASIS* basis, const int j);

  /*!
    index of first (leftmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  IntervalIndex<IBASIS> first_wavelet(const IBASIS* basis, const int j);

  /*!
    index of last (rightmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  IntervalIndex<IBASIS> last_wavelet(const IBASIS* basis, const int j);

  /*!
    index of first function with type e
    (mainly for TensorProductBasis)
  */
  template <class IBASIS>
  IntervalIndex<IBASIS> first_index(const IBASIS* basis, const int j, const int e);

  /*!
    index of last function with type e
    (mainly for TensorProductBasis)
  */
  template <class IBASIS>
  IntervalIndex<IBASIS> last_index(const IBASIS* basis, const int j, const int e);

  /*!
    number of first (leftmost) generator on level j0
  */
  template <class IBASIS>
  int first_generator_num(const IBASIS* basis);

  /*!
    number of last (rightmost) generator on level j0
  */
  template <class IBASIS>
  int last_generator_num(const IBASIS* basis);

  /*!
    number of first (leftmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  int first_wavelet_num(const IBASIS* basis, const int j);

  /*!
    number of last (rightmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  int last_wavelet_num(const IBASIS* basis, const int j);

}

// include implementation
#include <interval/i_index.cpp>

#endif

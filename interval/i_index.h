// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
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
    
    //! assignment
    IntervalIndex& operator = (const IntervalIndex& lambda);

    //! check equality
    bool operator == (const IntervalIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const IntervalIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    IntervalIndex& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const IntervalIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const IntervalIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    int j() const { return j_; }

    //! type e
    int e() const { return e_; }

    //! translation index k
    int k() const { return k_; }

    //! underlying basis
    const IBASIS* basis() const { return basis_; }

  protected:
    //! scale, type, translation
    int j_, e_, k_;

    //! pointer to the underlying interval basis
    const IBASIS* basis_;
  };

  //! stream output
  template <class IBASIS>
  inline std::ostream& operator << (std::ostream& os, const IntervalIndex<IBASIS>& lambda)
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
}

// include implementation
#include <interval/i_index.cpp>

#endif

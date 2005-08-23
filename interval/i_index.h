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
    A (template) wavelet index class for bases on [0,1] like those from [DKU],[DS].
    We require the template parameter class IBASIS to provide the routines
    j0(), DeltaLmin() and DeltaRmax(j).
  */
  template <class IBASIS>
  class IIndex
  {
  public:
    /*!
      constructor with given interval basis
      (also serves as a default constructor, but yields an invalid index then
      because the underlying interval basis must be specified to work correctly)
     */
    IIndex();

    //! copy constructor
    IIndex(const IIndex& lambda);
  
    //! constructor with specified parameters
    IIndex(const int j, const int e, const int k);
    
    //! assignment
    IIndex& operator = (const IIndex& lambda);

    //! check equality
    bool operator == (const IIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const IIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    IIndex& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const IIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const IIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    int j() const { return j_; }

    //! type e
    int e() const { return e_; }

    //! translation index k
    int k() const { return k_; }

  protected:
    //! scale, type, translation
    int j_, e_, k_;
  };

  //! stream output
  template <class IBASIS>
  inline std::ostream& operator << (std::ostream& os, const IIndex<IBASIS>& lambda)
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
  IIndex<IBASIS> first_generator(const int j);

  /*!
    index of last (rightmost) generator on level j >= j0
  */
  template <class IBASIS>
  IIndex<IBASIS> last_generator(const int j);

  /*!
    index of first (leftmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  IIndex<IBASIS> first_wavelet(const int j);

  /*!
    index of last (rightmost) wavelet on level j >= j0
  */
  template <class IBASIS>
  IIndex<IBASIS> last_wavelet(const int j);
}

// include implementation
#include <interval/i_index.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_DKU_INDEX_H
#define _WAVELETTL_DKU_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  template <int d, int dT>
  class DKUBasis;

  //! wavelet index class for bases on [0,1] from [DKU]
  template <int d, int dT>
  class DKUIndex
  {
  public:
    /*!
      constructor with given interval basis
      (also serves as a default constructor, but yields an invalid index then
      because the underlying interval basis must be specified to work correctly)
     */
    DKUIndex(const DKUBasis<d, dT>* basis = 0);

    //! copy constructor
    DKUIndex(const DKUIndex& lambda);
  
    //! constructor with specified parameters
    DKUIndex(const int j, const int e, const int k, const DKUBasis<d, dT>* basis);

    //! assignment
    DKUIndex& operator = (const DKUIndex& lambda);

    //! check equality
    bool operator == (const DKUIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const DKUIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    DKUIndex& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const DKUIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const DKUIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    int j() const { return j_; }

    //! type e
    int e() const { return e_; }

    //! translation index k
    int k() const { return k_; }
    
    //! underlying basis
    const DKUBasis<d, dT>* basis() const { return basis_; }

  protected:
    //! scale, type, translation
    int j_, e_, k_;

    //! pointer to the underlying DKU basis
    const DKUBasis<d, dT>* basis_;
  };

  //! stream output
  template <int d, int dT>
  inline std::ostream& operator << (std::ostream& os, const DKUIndex<d, dT>& lambda)
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
}

// include implementation
#include <interval/dku_index.cpp>

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_R_MW_INDEX_H
#define _WAVELETTL_R_MW_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  //! wavelet index class for (2-)multiwavelet bases over \mathbb R
  class RMWIndex
  {
  public:
    //! default constructor: yields unscaled generator index
    RMWIndex() : j_(0), e_(0), c_(0), k_(0) {}

    //! copy constructor
    RMWIndex(const RMWIndex& lambda);
  
    //! constructor with specified parameters
    RMWIndex(const int j, const int e, const int c, const int k);

    //! assignment
    RMWIndex& operator = (const RMWIndex& lambda);

    //! check equality
    bool operator == (const RMWIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const RMWIndex& lambda) const
    { return !(*this == lambda); }

    //! lexicographic order <
    bool operator < (const RMWIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const RMWIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    int j() const { return j_; }

    //! type e
    int e() const { return e_; }

    //! component c
    int c() const { return c_; }

    //! translation index k
    int k() const { return k_; }
  
  protected:
    //! scale, type, component, translation
    int j_, e_, c_, k_;
  };

  //! stream output
  inline std::ostream& operator << (std::ostream& os, const RMWIndex& lambda)
  {
    using namespace std;
    os << "("
       << lambda.j()
       << ","
       << lambda.e()
       << ","
       << lambda.c()
       << ","
       << lambda.k()
       << ")";
    return os;
  }
}

// include implementation
#include <Rd/r_mw_index.cpp>

#endif

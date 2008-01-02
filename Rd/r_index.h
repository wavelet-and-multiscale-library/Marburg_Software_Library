// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_R_INDEX_H
#define _WAVELETTL_R_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  //! wavelet index class for bases over \mathbb R
  class RIndex
  {
  public:
    //! default constructor: yields unscaled generator index
    RIndex() : j_(0), e_(0), k_(0) {}

    //! copy constructor
    RIndex(const RIndex& lambda);
  
    //! constructor with specified parameters
    RIndex(const int j, const int e, const int k);

    //! assignment
    RIndex& operator = (const RIndex& lambda);

    //! check equality
    bool operator == (const RIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const RIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    /*!
      preincrement; since the range for k is unbounded, we just
      have to increment k here
     */
    RIndex& operator ++ ()
    { ++k_; return *this; }
    
    //! lexicographic order <
    bool operator < (const RIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const RIndex& lambda) const
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
  inline std::ostream& operator << (std::ostream& os, const RIndex& lambda)
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
#include <Rd/r_index.cpp>

#endif

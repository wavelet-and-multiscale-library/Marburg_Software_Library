// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_R_Q_INDEX_H
#define _WAVELETTL_R_Q_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  //! index class for quarklet frames over \mathbb R
  class RQIndex
  {
  public:
    //! default constructor: yields unscaled generator index
    RQIndex() : p_(0), j_(0), e_(0), k_(0) {}

    //! copy constructor
    RQIndex(const RQIndex& lambda);
  
    //! constructor with specified parameters
    RQIndex(const int p, const int j, const int e, const int k);

    //! assignment
    RQIndex& operator = (const RQIndex& lambda);

    //! check equality
    bool operator == (const RQIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const RQIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    /*!
      preincrement; since the range for k is unbounded, we just
      have to increment k here
     */
    RQIndex& operator ++ ()
    { ++k_; return *this; }
    
    //! lexicographic order <
    bool operator < (const RQIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const RQIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! polynomial degree p 
    int p() const { return p_; }

    //! scale j
    int j() const { return j_; }

    //! type e
    int e() const { return e_; }

    //! type index type
    typedef int type_type;

    //! translation index k
    int k() const { return k_; }
  
    //! translation index type
    typedef int translation_type;

  protected:
    //! polynomial degree, scale, type, translation
    int p_, j_, e_, k_;
  };

  //! stream output
  inline std::ostream& operator << (std::ostream& os, const RQIndex& lambda)
  {
    using namespace std;
    os << "("
       << lambda.p()
       << ","
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
#include <Rd/r_q_index.cpp>

#endif

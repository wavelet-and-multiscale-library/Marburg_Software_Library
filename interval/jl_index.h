// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_JL_INDEX_H
#define _WAVELETTL_JL_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

namespace WaveletTL
{
  /*!
    Index class for the multiwavelet basis in [0,1] from [JL],
    modelling a multiindex of the form lambda=(j,e,c,k).

    We implicitly assume that j0==1.
  */
  class JLIndex
  {
  public:
    //! default constructor
    JLIndex();

    //! copy constructor
    JLIndex(const JLIndex& lambda);

    //! copy index from const pointer
    JLIndex(const JLIndex* lambda);

    //! type index type
    typedef int type_type;

    //! component type
    typedef int component_type;

    //! translation index type
    typedef int translation_type;

    //! constructor with specified parameters
    JLIndex(const int j, const type_type e, const component_type c, const translation_type k);

    //! assignment
    JLIndex& operator = (const JLIndex& lambda);

    //! check equality
    bool operator == (const JLIndex& lambda) const;

    //! check non-equality
    bool operator != (const JLIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    JLIndex& operator ++ ();
    
    //! lexicographic order <
    bool operator < (const JLIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const JLIndex& lambda) const
    { return (*this < lambda || *this == lambda); }
    
    //! scale j
    const int j() const { return j_; }
    
    //! type e
    const type_type& e() const { return e_; }
    
    //! component type c
    const component_type& c() const { return c_; }

    //! translation index k
    const translation_type& k() const { return k_; }

  protected:  
    
    //! scale, type, component, translation
    int j_, e_, c_, k_;
    
  };

  //! stream output
  std::ostream& operator << (std::ostream& os, const JLIndex& lambda)
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

  /*!
    index of first (leftmost) generator on level j >= j0
  */
  JLIndex first_generator(const int j);

  /*!
    index of last (rightmost) generator on level j >= j0
  */
  JLIndex last_generator(const int j);

  /*!
    index of first (leftmost) wavelet on level j >= j0
  */
  JLIndex first_wavelet(const int j);

  /*!
    index of last (rightmost) wavelet on level j >= j0
  */
  JLIndex last_wavelet(const int j);
}

// include implementation
#include <interval/jl_index.cpp>

#endif

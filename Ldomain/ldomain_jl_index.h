// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_JL_INDEX_H
#define _WAVELETTL_LDOMAIN_JL_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

#include <utils/multiindex.h>

using MathTL::MultiIndex;

namespace WaveletTL
{
  class LDomainJLBasis;
  
  /*!
    an index class for LDomainJLBasis;
    we implicitly assume that j0==1.
  */
  class LDomainJLIndex
  {
  public:
    //! type index type
    typedef MultiIndex<int,2> type_type;
    
    //! component index type
    typedef MultiIndex<int,2> component_type;
    
    //! translation index type
    typedef MultiIndex<int,2> translation_type;
    
    //! default constructor
    LDomainJLIndex();
    
    //! copy constructor
    LDomainJLIndex(const LDomainJLIndex& lambda);
    
    //! constructor with given j,e,c,k
    LDomainJLIndex(const int j,
		   const type_type& e,
		   const component_type& c,
		   const translation_type& k);
    
    //! assignment
    LDomainJLIndex& operator = (const LDomainJLIndex& lambda);

    //! check equality
    bool operator == (const LDomainJLIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const LDomainJLIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    LDomainJLIndex& operator ++ ();

    //! lexicographic order <
    bool operator < (const LDomainJLIndex& lambda) const;
    
    //! lexicographic order <=
    bool operator <= (const LDomainJLIndex& lambda) const
    { return (*this < lambda || *this == lambda); }
    
    //! scale j
    const int j() const { return j_; }
    
    //! type e
    const type_type& e() const { return e_; }
    
    //! component c
    const component_type& c() const { return c_; }
    
    //! translation index k
    const translation_type& k() const { return k_; }
    
  protected:
    //! scale
    int j_;
    
    //! type
    type_type e_;
    
    //! component
    component_type c_;
    
    //! translation
    translation_type k_;
  };
  
  //! stream output
  std::ostream& operator << (std::ostream& os,
			     const LDomainJLIndex& lambda)
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
  
  //! index of first generator on level j >= j0
  LDomainJLIndex
  first_generator(const int j);

  //! index of last generator on level j >= j0
  LDomainJLIndex
  last_generator(const int j);

  //! index of first wavelet on level j >= j0 
  LDomainJLIndex
  first_wavelet(const int j);
  
  //! index of first wavelet with a given type on level j >= j0
  LDomainJLIndex
  first_wavelet(const int j,
 		const LDomainJLIndex::type_type& e);
  
  //! index of last wavelet on level j >= j0
  LDomainJLIndex
  last_wavelet(const int j);

  //! quickly decide whether a given index would be valid
  bool index_is_valid(const int j, const int e0, const int e1,
		      const int c0, const int c1, const int k0, const int k1);

}

#include <Ldomain/ldomain_jl_index.cpp>

#endif

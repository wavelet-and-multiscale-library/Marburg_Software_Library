// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2006                                            |
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
    
    //! constructor with given j,e,c,k
    LDomainJLIndex(const int j,
		   const type_type& e,
		   const component_type& c,
		   const translation_type& k);
    
    //! check equality
    bool operator == (const LDomainJLIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const LDomainJLIndex& lambda) const
    { return !(*this == lambda); }

    //! preincrement
    LDomainJLIndex& operator ++ ();

/*     //! lexicographic order < */
/*     bool operator < (const LDomainIndex& lambda) const; */

/*     //! lexicographic order <= */
/*     bool operator <= (const LDomainIndex& lambda) const */
/*     { return (*this < lambda || *this == lambda); } */
    
    //! scale j
    const int j() const { return j_; }
    
    //! type e
    const type_type& e() const { return e_; }
    
    //! component c
    const component_type& c() const { return c_; }
    
    //! translation index k
    const translation_type& k() const { return k_; }
    
/*     /\*! */
/*       By construction, the overall wavelet index set is ordered, so that */
/*       there exists a bijective mapping into the positive integers. */
/*       This routine returns the "number" of the current index, starting with 0 */
/*       for the first generator on the coarsest level. If the current index is a */
/*       generator on a higher level j, the number 0 corresponds to the first generator */
/*       on the level j. */
/*     *\/ */
/*     const int number() const; */
    
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

/*   /\*! */
/*     index of first wavelet on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS> */
/*   LDomainIndex<IBASIS> */
/*   first_wavelet(const LDomainBasis<IBASIS>* basis, const int j); */
  
/*   /\*! */
/*     index of first wavelet with a given type on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS> */
/*   LDomainIndex<IBASIS> */
/*   first_wavelet(const LDomainBasis<IBASIS>* basis, */
/* 		const int j, */
/* 		const typename LDomainIndex<IBASIS>::type_type& e); */
  
  //! index of last wavelet on level j >= j0
  LDomainJLIndex
  last_wavelet(const int j);
}

#include <Ldomain/ldomain_jl_index.cpp>

#endif

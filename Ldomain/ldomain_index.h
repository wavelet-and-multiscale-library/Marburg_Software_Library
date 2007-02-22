// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_LDOMAIN_INDEX_H
#define _WAVELETTL_LDOMAIN_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

#include <utils/multiindex.h>

using MathTL::MultiIndex;

namespace WaveletTL
{
  template <class IBASIS> class LDomainBasis;
  
  /*!
    An index class for composite wavelet bases of the L-shaped domain in R^2.
    In fact, this is essentially a lexicographically ordered multiindex (j,e,p,k), with
    * scale j
    * type e ((0,0),(0,1),(1,0) or (1,1))
    * (logical) patch p (0<=p<=4)
    * translation index k

    The logical patches are:
      p=0 (-1,0)x( 0,1)
      p=1 (-1,0)x(-1,0)
      p=2 ( 0,1)x(-1,0)
      p=3 (-1,0)x  {0}
      p=4   {0} x(-1,0)
  */
  template <class IBASIS>
  class LDomainIndex
  {
  public:
    //! type index type
    typedef MultiIndex<int,2> type_type;
    
    //! translation index type
    typedef MultiIndex<int,2> translation_type;

    /*!
      constructor with a given L-domain basis
      (also serves as a default constructor, but yields an invalid index
      in this case, because the underlying bases must be specified to work correctly)
    */
    LDomainIndex(const LDomainBasis<IBASIS>* basis = 0);

    //! constructor with given j,e,p,k
    LDomainIndex(const int j,
		 const type_type& e,
		 const int p,
		 const translation_type& k,
		 const LDomainBasis<IBASIS>* basis);

    //! copy constructor
    LDomainIndex(const LDomainIndex& lambda);

    //! copy index from const pointer
    LDomainIndex(const LDomainIndex* lambda);

  
    //! check equality
    bool operator == (const LDomainIndex& lambda) const;
    
    //! check non-equality
    inline bool operator != (const LDomainIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    LDomainIndex& operator ++ ();

    //! lexicographic order <
    bool operator < (const LDomainIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const LDomainIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    const int j() const { return j_; }

    //! type e
    const type_type& e() const { return e_; }

    //! patch p
    const int p() const { return p_; }

    //! translation index k
    const translation_type& k() const { return k_; }

    /*!
      By construction, the overall wavelet index set is ordered, so that
      there exists a bijective mapping into the positive integers.
      This routine returns the "number" of the current index, starting with 0
      for the first generator on the coarsest level. If the current index is a
      generator on a higher level j, the number 0 corresponds to the first generator
      on the level j.
    */
    const int number() const;
    
  protected:
    //! pointer to the underlying basis
    const LDomainBasis<IBASIS>* basis_;

    //! scale
    int j_;
    
    //! type
    MultiIndex<int,2> e_;

    //! patch
    int p_;
    
    //! translation
    MultiIndex<int,2> k_;

  };

  //! stream output
  template <class IBASIS>
  inline std::ostream& operator << (std::ostream& os,
				    const LDomainIndex<IBASIS>& lambda)
  {
    using namespace std;
    os << "("
       << lambda.j()
       << ","
       << lambda.e()
       << ","
       << lambda.p()
       << ","
       << lambda.k()
       << ")";
    return os;
  }

  /*!
    index of first generator on level j >= j0
  */
  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_generator(const LDomainBasis<IBASIS>* basis, const int j);

  /*!
    index of last generator on level j >= j0
  */
  template <class IBASIS>
  LDomainIndex<IBASIS>
  last_generator(const LDomainBasis<IBASIS>* basis, const int j);

  /*!
    index of first wavelet on level j >= j0
  */
  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_wavelet(const LDomainBasis<IBASIS>* basis, const int j);
  
  /*!
    index of first wavelet with a given type on level j >= j0
  */
  template <class IBASIS>
  LDomainIndex<IBASIS>
  first_wavelet(const LDomainBasis<IBASIS>* basis,
		const int j,
		const typename LDomainIndex<IBASIS>::type_type& e);
  
  /*!
    index of last wavelet on level j >= j0
  */
  template <class IBASIS>
  LDomainIndex<IBASIS>
  last_wavelet(const LDomainBasis<IBASIS>* basis, const int j);
}

#include <Ldomain/ldomain_index.cpp>

#endif

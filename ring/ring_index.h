// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_RING_INDEX_H
#define _WAVELETTL_RING_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

#include <utils/multiindex.h>

using MathTL::MultiIndex;

namespace WaveletTL
{
  template <int d, int dt, int s0, int s1> class RingBasis;
  
  /*!
    An index class for tensor product wavelet bases on a
    ring-shaped domain in R^2.
  */
  template <int d, int dt, int s0, int s1>
  class RingIndex
  {
  public:
    //! type index type
    typedef MultiIndex<int,2> type_type;
    
    //! translation index type
    typedef MultiIndex<int,2> translation_type;

    //! default constructor
    RingIndex();

    /*!
      constructor with given j,e,k
    */
    RingIndex(const int j,
	      const type_type& e,
	      const translation_type& k);

    //! copy constructor
    RingIndex(const RingIndex& lambda);

    //! copy index from const pointer
    RingIndex(const RingIndex* lambda);

    //! assignment
    RingIndex& operator = (const RingIndex& lambda);

    //! check equality
    bool operator == (const RingIndex& lambda) const;
    
    //! check non-equality
    inline bool operator != (const RingIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    RingIndex& operator ++ ();

    //! lexicographic order <
    bool operator < (const RingIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const RingIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    const int j() const { return j_; }

    //! type e
    const type_type& e() const { return e_; }

    //! translation index k
    const translation_type& k() const { return k_; }

  protected:
    //! scale
    int j_;
    
    //! type
    type_type e_;

    //! translation
    translation_type k_;

  };

  //! stream output
  template <int d, int dt, int s0, int s1>
  inline std::ostream& operator << (std::ostream& os,
				    const RingIndex<d,dt,s0,s1>& lambda)
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

/*   /\*! */
/*     index of first generator on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   CubeIndex<IBASIS,DIM,CUBEBASIS> */
/*   first_generator(const CUBEBASIS* basis, const int j); */

/*   /\*! */
/*     index of last generator on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   CubeIndex<IBASIS,DIM,CUBEBASIS> */
/*   last_generator(const CUBEBASIS* basis, const int j); */

/*   /\*! */
/*     index of first wavelet on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   CubeIndex<IBASIS,DIM,CUBEBASIS> */
/*   first_wavelet(const CUBEBASIS* basis, const int j); */
  
/*   /\*! */
/*     index of last wavelet on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   CubeIndex<IBASIS,DIM,CUBEBASIS> */
/*   last_wavelet(const CUBEBASIS* basis, const int j); */

/*   /\*! */
/*     number of first generator on level j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   const int */
/*   first_generator_num(const CUBEBASIS* basis); */

/*   /\*! */
/*     number of last generator on level j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   const int */
/*   last_generator_num(const CUBEBASIS* basis); */

/*   /\*! */
/*     number of first wavelet on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   const int */
/*   first_wavelet_num(const CUBEBASIS* basis, const int j); */
  
/*   /\*! */
/*     number of last wavelet on level j >= j0 */
/*   *\/ */
/*   template <class IBASIS, unsigned int DIM, class CUBEBASIS> */
/*   const int */
/*   last_wavelet_num(const CUBEBASIS* basis, const int j); */

}

#include <ring/ring_index.cpp>

#endif

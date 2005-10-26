// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_CUBE_INDEX_H
#define _WAVELETTL_CUBE_INDEX_H

#include <iostream>
using std::cout;
using std::endl;

#include <utils/multiindex.h>

using MathTL::MultiIndex;

namespace WaveletTL
{
  template <class IBASIS, unsigned int DIM> class CubeBasis;
  
  /*!
    An index class for tensor product wavelet bases over the d-dimensional
    unit cube [0,1]^d (or mapped versions thereof).
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS = CubeBasis<IBASIS,DIM> >
  class CubeIndex
  {
  public:
    //! type index type
    typedef MultiIndex<unsigned int,DIM> type_type;
    
    //! translation index type
    typedef MultiIndex<int,DIM> translation_type;

    /*!
      constructor with a given cube basis
      (also serves as a default constructor, but yields an invalid index pair
      in this case, because the underlying bases must be specified to work correctly)
    */
    CubeIndex(const CUBEBASIS* basis = 0);

    /*!
      constructor with given j,e,k
    */
    CubeIndex(const int j,
	      const type_type& e,
	      const translation_type& k,
	      const CUBEBASIS* basis);

    //! copy constructor
    CubeIndex(const CubeIndex& lambda);
  
    //! check equality
    bool operator == (const CubeIndex& lambda) const;
    
    //! check non-equality
    inline bool operator != (const CubeIndex& lambda) const
    { return !(*this == lambda); }
    
    //! preincrement
    CubeIndex& operator ++ ();

    //! lexicographic order <
    bool operator < (const CubeIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const CubeIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! scale j
    const int j() const { return j_; }

    //! type e
    const type_type& e() const { return e_; }

    //! translation index k
    const translation_type& k() const { return k_; }

  protected:
    //! pointer to the underlying basis
    const CUBEBASIS* basis_;

    //! scale
    int j_;
    
    //! type
    MultiIndex<unsigned int,DIM> e_;

    //! translation
    MultiIndex<int,DIM> k_;
  };

  //! stream output
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  inline std::ostream& operator << (std::ostream& os,
				    const CubeIndex<IBASIS,DIM,CUBEBASIS>& lambda)
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
    index of first generator on level j >= j0
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  first_generator(const CUBEBASIS* basis, const int j);

  /*!
    index of last generator on level j >= j0
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  last_generator(const CUBEBASIS* basis, const int j);

  /*!
    index of first wavelet on level j >= j0
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  first_wavelet(const CUBEBASIS* basis, const int j);
  
  /*!
    index of last wavelet on level j >= j0
  */
  template <class IBASIS, unsigned int DIM, class CUBEBASIS>
  CubeIndex<IBASIS,DIM,CUBEBASIS>
  last_wavelet(const CUBEBASIS* basis, const int j);
}

#include <cube/cube_index.cpp>

#endif

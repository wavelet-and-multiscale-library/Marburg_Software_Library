// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MULTIINDEX_H
#define _MATHTL_MULTIINDEX_H

#include <iostream>
#include <set>
#include <utils/fixed_array1d.h>

namespace MathTL
{
  /*!
    Template class for (homogeneous) multiindices of a priori known length.
    The class is designed to be used as index class in a std::map<I,.>.
  */
  template <class I, unsigned int DIMENSION>
  class MultiIndex
    : public FixedArray1D<I, DIMENSION>
  {
  public:
    /*!
      default constructor, yields a zero multiindex
    */
    MultiIndex();

    /*!
      copy constructor
    */
    MultiIndex(const MultiIndex& lambda);

    /*!
      constructor from a single index, this is only allowed for DIMENSION == 1
    */
    explicit MultiIndex(const I& i0);

    /*!
      constructor from two single indices, this is only allowed for DIMENSION == 2
    */
    MultiIndex(const I& i0, const I& i1);
  
    //! check equality
    bool operator == (const MultiIndex& lambda) const;

    //! check non-equality
    inline bool operator != (const MultiIndex& lambda) const
    { return !(*this == lambda); }

    //! lexicographic order <
    bool operator < (const MultiIndex& lambda) const;

    //! lexicographic order <=
    bool operator <= (const MultiIndex& lambda) const
    { return (*this < lambda || *this == lambda); }
    
  };

  /*!
    For two multiindices \alpha and \beta, return the set
    of all multiindices \gamma such that the componentwise property
      \alpha \le \gamma \le \beta
    holds.
  */
  template<class I, unsigned int DIMENSION>
  std::set<MultiIndex<I, DIMENSION> >
  cuboid_indices(const MultiIndex<I, DIMENSION>& alpha,
		 const MultiIndex<I, DIMENSION>& beta);

  /*!
    Compute all multiindices \alpha\in\mathbb N^d with degree |\alpha|=k
  */
  template <unsigned int DIMENSION>
  std::set<MultiIndex<unsigned int, DIMENSION> >
  degree_indices(const unsigned int k);
  
  //! stream output
  template<class I, unsigned int DIMENSION>
  inline std::ostream&
  operator << (std::ostream& os, const MultiIndex<I, DIMENSION>& lambda)
  {
    using namespace std;
    os << "(";
    for (unsigned int i(0); i < DIMENSION; i++)
      {
	os << lambda[i];
	if (i < DIMENSION-1)
	  os << ",";
	else
	  os << ")";
      }
    
    return os;
  }
}

#include <utils/multiindex.cpp>

#endif

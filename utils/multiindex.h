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
#include <utils/array1d.h>

namespace MathTL
{
  /*!
    Template class for (homogeneous) multiindices of a priori known length.
    The class is designed to be used as index class in a std::map<I,.>.
  */
  template <class I, unsigned int DIMENSION>
  class MultiIndex
    : public Array1D<I>
  {
  public:
    /*!
      default constructor, yields a zero multiindex
    */
    MultiIndex();

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

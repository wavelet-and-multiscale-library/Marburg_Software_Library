// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _WAVELETTL_PERIODIC_H
#define _WAVELETTL_PERIODIC_H

#include <iostream>

namespace WaveletTL
{
  /*!
    Template class for a periodic, biorthogonal wavelet basis on the unit interval [0,1].
    The template parameters provide the biorthogonal masks of a refinable function
    on the real line.
    If the second template parameter is omitted, we assume that the given refinable
    function induces an orthonormal wavelet basis.

    A univariate mask should have at least the signature of a LaurentPolynomial<double>,
    i.e., have the usual iterator classes.
  */
  template <class PRIMALMASK, class DUALMASK = PRIMALMASK>
    class PeriodicBasis
    {
      public:
      /*!
	Wavelet index class for periodic wavelet bases on [0,1]
       */
      class Index
      {
	public:
	//! default constructor: yields unscaled generator index
	Index() : j_(0), e_(0), k_(0) {}
	
	//! copy constructor
	Index(const Index& lambda);
	
	//! constructor with specified parameters
	Index(const int j, const int e, const int k);
	
	//! assignment
	Index& operator = (const Index& lambda);

	//! check equality
	bool operator == (const Index& lambda) const;

	//! check non-equality
	inline bool operator != (const Index& lambda) const
	{ return !(*this == lambda); }
	
	//! preincrement
	Index& operator ++ ();

	//! scale j
	inline int j() const { return j_; }
	
	//! type e (e=0: generator, e=1: wavelet)
	inline int e() const { return e_; }
	
	//! translation index k
	inline int k() const { return k_; }

	//! lexicographic order <
	bool operator < (const Index& lambda) const;
	
	//! lexicographic order <=
	bool operator <= (const Index& lambda) const
	{ return (*this < lambda || *this == lambda); }

	protected:
	//! scale, type, geometric location
	int j_, e_, k_;
      };
    };
}

// include implementation
#include <interval/periodic.cpp>

#endif

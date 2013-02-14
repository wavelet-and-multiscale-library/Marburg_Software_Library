// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Ulrich Friedrich                   |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MULTIINDEX_H
#define _MATHTL_MULTIINDEX_H

#include <iostream>
#include <set>
#include <map>
#include <utils/fixed_array1d.h>

namespace MathTL
{
  /*!
   * Template class for (homogeneous) multiindices of a priori known length with values in \N^d.
   * The class is designed to be used as index class in a std::map<I,.>.
   * Values in \Z^d are not prohibited, but also not supported. 
   * The ordering and numbering will be wrong in this case!
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

    /*
     * MultiIndices are ordered same as \N^dim is ordered.
     * This is not the lexicographic order!
     * This does not work for multiindices with negative entries.
     * Explicitly: indices are ordered primarily by their 1 norm, secondly lexicographically.
     */
    bool operator < (const MultiIndex& lambda) const;

    /*
     * Lexicographical ordering
     */
    bool lex (const MultiIndex& lambda) const;

    //! ordered by distance from 0
    bool operator <= (const MultiIndex& lambda) const
    { return (*this < lambda || *this == lambda); }

    //! ordered by distance from 0
    bool operator >= (const MultiIndex& lambda) const
    { return (lambda < *this || *this == lambda); }
 
    /*
     * Preincrement. Works only for nonnegative entries!
     * Ordering given by numbering of \N^dim, i.e.
     * (0,0,0),(0,0,1),(0,1,0),(1,0,0),(0,0,2),(0,1,1),(0,2,0),(1,0,1),(1,1,0),(2,0,0),(0,0,3),...
     * This ordering is used for the operator <.
     */
    MultiIndex& operator ++ ();
    
	// Postincrement. Works only for nonnegative entries!
    MultiIndex operator ++ (int)    
    {
    	MultiIndex temp (*this);
    	++*this;
    	return temp;
    };

    /*
     * Return Number of MulitiIndex with respect to the ordering "<"
     * Counting begins at 0.
     */
    unsigned long int number();
    
  };


  /*
   * Returns the mapping from the multiindices to the natural numbers according 
   * to the ordering <.
   * Output contains multiindices ranging from 0 to lambda, excluding lambda.
   */
  template <class I, unsigned int DIMENSION>
  std::map<MultiIndex<I, DIMENSION>,int>
  indexmapping(const MultiIndex<I, DIMENSION>& lambda);

  /*!
    For two multiindices \alpha and \beta, return the set
    of all multiindices \gamma such that the componentwise property
      \alpha \le \gamma \le \beta
    holds.
  */
  template <class I, unsigned int DIMENSION>
  std::set<MultiIndex<I, DIMENSION> >
  cuboid_indices(const MultiIndex<I, DIMENSION>& alpha,
		 const MultiIndex<I, DIMENSION>& beta);

  /*!
   * degree of a multiindex \alpha\in\mathbb N^d, i.e.,
   * the \ell_1 norm 
   * (I should be "int" or "unsigned int")
  */
  template <class I, unsigned int DIMENSION>
  unsigned int multi_degree(const MultiIndex<I, DIMENSION>& alpha);

  /*!
   * factorial of a multiindex \alpha\in\mathbb N^d, i.e.,
   * alpha ! = alpha_1! * ... * alpha_n!
   * (I should be "int" or "unsigned int")
  */
  template <class I, unsigned int DIMENSION>
  unsigned int multi_factorial(const MultiIndex<I, DIMENSION>& alpha);

  /*!
   * \beta-th power of a multiindex \alpha, i.e.,
   * alpha^beta = alpha_1^beta_1 * ... * alpha_n^beta_n
   * (I should be "int" or "unsigned int")
  */
  template <class I, unsigned int DIMENSION>
  int multi_power(const MultiIndex<I, DIMENSION>& alpha,
		  const MultiIndex<I, DIMENSION>& beta);
  
  /*!
   * binomial \beta over \alpha of two multiindices, i.e.,
   * (alpha beta) = (alpha_1 beta_1) * ... * (alpha_n beta_n)
  */
  template <class I, unsigned int DIMENSION>
  int multi_binomial(const MultiIndex<I, DIMENSION>& beta,
		     const MultiIndex<I, DIMENSION>& alpha);

  /*!
    Compute all multiindices \alpha\in\mathbb N^d with degree |\alpha|=k
    (I should be "int" or "unsigned int")
  */
  template <class I, unsigned int DIMENSION>
  std::set<MultiIndex<I, DIMENSION> >
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

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_DIFFERENCES_H
#define _MATHTL_DIFFERENCES_H

#include <iostream>
#include <algebra/infinite_vector.h>

namespace MathTL
{
  /*!
    Evaluate a k-th forward difference on an infinite vector of
    function values
    (which we assume to correspond to a uniformly sampled
    mapping on the integers).
  */
  template <unsigned int K>
  inline
  InfiniteVector<double, int> forward_difference(const InfiniteVector<double, int>& a)
  {
    return forward_difference<1>(forward_difference<K-1>(a));
  }

  /*!
    the 0-th forward difference is the identity
   */
  template <>
  inline
  InfiniteVector<double, int> forward_difference<0>(const InfiniteVector<double, int>& a)
  {
    return InfiniteVector<double, int>(a);
  }

  /*!
    forward difference
  */
  template <>
  InfiniteVector<double, int> forward_difference<1>(const InfiniteVector<double, int>& a)
  {
    InfiniteVector<double, int> r;

    if (!a.empty())
      {
 	r[a.begin().index()-1] = *(a.begin());
	
	for (int k = a.begin().index(); k <= a.rbegin().index(); k++)
	  {
	    double help = a[k+1]-a[k];
	    if (help != 0)
	      r[k] = help;
	  }
      }
	
    return r;
  }

  /*!
    Evaluate a k-th backward difference on an infinite vector of function values.
  */
  template <unsigned int K>
  inline
  InfiniteVector<double, int> backward_difference(const InfiniteVector<double, int>& a)
  {
    return backward_difference<1>(backward_difference<K-1>(a));
  }

  /*!
    the 0-th backward difference is the identity
  */
  template <>
  inline
  InfiniteVector<double, int> backward_difference<0>(const InfiniteVector<double, int>& a)
  {
    return InfiniteVector<double, int>(a);
  }

  /*!
    first backward difference
  */
  template <>
  InfiniteVector<double, int> backward_difference<1>(const InfiniteVector<double, int>& a)
  {
    InfiniteVector<double, int> r;

    if (!a.empty())
      {
 	r[a.rbegin().index()+1] = -*(a.rbegin());
	
	for (int k = a.begin().index(); k <= a.rbegin().index(); k++)
	  {
	    double help = a[k]-a[k-1];
	    if (help != 0)
	      r[k] = help;
	  }
      }
	
    return r;
  }
}

#endif

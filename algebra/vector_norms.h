// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_VECTOR_NORMS_H
#define _MATHTL_VECTOR_NORMS_H

#include <algorithm>
#include <cmath>
#include <iostream>

namespace MathTL
{
  /*
    computation of different (quasi-)norms for generic VECTOR classes
    with a standard signature like std::vector<T>

    You can use template specialization to speed up the norm computation
    for specific VECTOR classes.
  */

  /*!
    maximum norm
   */
  template <class VECTOR>
  const double linfty_norm(const VECTOR& v)
  {
    double r(0.0);

    for (typename VECTOR::const_iterator it(v.begin()), itend(v.end());
	 it != itend; ++it)
      r = std::max(fabs(*it), r);

    return r;
  }

  /*!
    l_1 norm
  */
  template <class VECTOR>
  const double l1_norm(const VECTOR& v)
  {
    double r(0.0);
    
    for (typename VECTOR::const_iterator it(v.begin()), itend(v.end());
	 it != itend; ++it)
      r += fabs(*it);

    return r;
  }

  /*!
    square of l_2 norm
  */
  template <class VECTOR>
  const double l2_norm_sqr(const VECTOR& v)
  {
    double r(0.0);

    for (typename VECTOR::const_iterator it(v.begin()), itend(v.end());
	 it != itend; ++it)
      r += *it * *it;

    return r;
  }

  /*!
    l_2 norm
  */
  template <class VECTOR>
  const double l2_norm(const VECTOR& v)
  {
    return sqrt(l2_norm_sqr(v));
  }

  /*!
    l_p (quasi-) norm
  */
  template <class VECTOR>
  const double lp_norm(const VECTOR& v, const double p)
  {
    double r(0.0);

    for (typename VECTOR::const_iterator it(v.begin()), itend(v.end());
	 it != itend; ++it)
      r += std::pow(fabs(*it), p);
    
    return std::pow(r, 1.0/p);
  }
}

#endif

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_VECTOR_NORMS_H
#define _MATHTL_VECTOR_NORMS_H

#include <algorithm>
#include <iostream>

namespace MathTL
{
  /*
    computation of different (quasi-)norms for generic VECTOR classes
    with a standard signatur like std::vector<T>

    We use the size() routine for determining the dimension of the vector.

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

    for (unsigned int i(0); i < v.size(); i++)
      r = std::max(fabs(v[i]), r);

    return r;
  }

  /*!
    l_1 norm
  */
  template <class VECTOR>
  const double l1_norm(const VECTOR& v)
  {
    double r(0.0);
    
    for (unsigned int i(0); i < v.size(); i++)
      r += fabs(v[i]);

    return r;
  }

  /*!
    square of l_2 norm
  */
  template <class VECTOR>
  const double norm_sqr(const VECTOR& v)
  {
    double r(0.0);

    for (unsigned int i(0); i < v.size(); i++)
      r += v[i] * v[i];

    return r;
  }

  /*!
    l_2 norm
  */
  template <class VECTOR>
  const double l2_norm(const VECTOR& v)
  {
    return sqrt(norm_sqr(v));
  }

  /*!
    l_p (quasi-) norm
  */
  template <class VECTOR>
  const double lp_norm(const VECTOR& v, const double p)
  {
    double r(0.0);

    for (unsigned int i(0); i < v.size(); i++)
      r += pow(fabs(v[i]), p);
    
    return pow(r, 1.0/p);
  }
}

#endif

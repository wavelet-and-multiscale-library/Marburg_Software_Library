// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_VECTOR_ARITHMETICS_H
#define _MATHTL_VECTOR_ARITHMETICS_H

#include <algorithm>
#include <cassert>

namespace MathTL
{
  /*
    diverse external arithmetics for generic VECTOR classes
    with a standard signatur like std::vector<T>

    We use the size() routine for determining the dimension of the vector.

    You can use template specialization to speed up computation
    for specific VECTOR classes.
  */

  /*!
    mean value of a vector
  */
  template <class VECTOR>
  double mean_value(const VECTOR& v)
  {
    assert(v.size() > 0);

    double r(0);

    typename VECTOR::const_iterator it(v.begin()), itend(v.end());
    while (it != itend)
      r += *it++;

    return r/v.size();
  }
}

#endif

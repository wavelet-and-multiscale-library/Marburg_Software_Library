// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_VECTOR_IO_H
#define _MATHTL_VECTOR_IO_H

#include <iostream>

namespace MathTL
{
  /*
    generic input/output routines for VECTOR classes with a standard
    signature like of std::vector<T>
    We use the size() routine for determining the dimension of the vector.
  */
  
  /*!
    generic Matlab-style stream output of the form
    [x_1 x_2 ... x_n]
  */
  template <class VECTOR>
  void print_vector(const VECTOR& v, std::ostream& os)
  {
    os << "[";
    if (v.size() > 0)
      {
	os << v[0];
	for (unsigned int i(1); i < v.size(); i++)
	  os << " " << v[i];
      }
    os << "]";
  }
}

#endif

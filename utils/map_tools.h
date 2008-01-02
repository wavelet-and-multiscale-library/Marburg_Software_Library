// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_MAP_TOOLS_H
#define _MATHTL_MAP_TOOLS_H

#include <cassert>
#include <map>

namespace MathTL
{
  /*!
    add two maps in O(N), result=factor1*m1+factor2*m2
  */
  template <class K, class C>
  void
  add_maps(const std::map<K,C>& m1, const std::map<K,C>& m2,
	   std::map<K,C>& result,
	   const double factor1 = 1.0, const double factor2 = 1.0);
  
  /*!
    stream output for an arbitrary map
  */
  template <class K, class C>
  std::ostream& operator << (std::ostream& os, const std::map<K,C>& m);
  
}

#include <utils/map_tools.cpp>

#endif

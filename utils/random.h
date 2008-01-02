// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2008                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_RANDOM_H
#define _MATHTL_RANDOM_H

#include <cmath>
#include <cstdlib>
#include <time.h>

namespace MathTL
{
  /*!
    generate uniformly distributed double random numbers in [a,b)
  */
  double random_double(const double a = 0, const double b = 1)
  {
    return a + (b-a) * (double(rand())/RAND_MAX);
  }

  /*!
    generate uniformly distributed integer random numbers in [a,b]
  */
  int random_integer(const int a = 0, const int b = 1)
  {
    return (int) (floor) (a + (b-a+1) * (double(rand())/RAND_MAX));
  } 
}

#endif

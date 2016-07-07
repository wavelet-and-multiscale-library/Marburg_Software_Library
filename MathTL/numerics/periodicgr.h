// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_PERIODICGR_H
#define _MATHTL_PERIODICGR_H

/*Base class for the Periodic Gramian equations u(t)=f(t) on [0,1] with periodic boundary 
 * conditions u(0)=u(1) and \int_0^1 u=0. 
 */
namespace MathTL
{
  class PeriodicGramianProblem
  {
  public:
    /*!
      virtual destructor
     */
    virtual ~PeriodicGramianProblem ();
    virtual double g(const double t) const = 0;
  };
}
#include <numerics/periodicgr.cpp>
#endif

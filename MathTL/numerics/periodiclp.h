// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner, Philipp Keding                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_PERIODICLP_H
#define _MATHTL_PERIODICLP_H

/*Base class for the problem -u''(t)=f(t) on [0,1] with periodic boundary 
 * conditions u(0)=u(1) and \int_0^1 u=0. 
 * (Can also be used for Periodic Gramian equations u(t)=f(t) on [0,1])@PHK*/
namespace MathTL
{
  class PeriodicLaplacianProblem
  {
  public:
    /*!
      virtual destructor
     */
    virtual ~PeriodicLaplacianProblem ();
    virtual double g(const double t) const = 0;
  };
}
#include <numerics/periodiclp.cpp>
#endif

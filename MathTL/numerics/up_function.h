// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_UP_FUNCTION_H
#define _MATHTL_UP_FUNCTION_H

#include <utils/function.h>
#include <utils/array1d.h>

namespace MathTL
{
  
  /*!
    It is well-known that the sequence of functions
    
      u_0(x) = \chi_{[-1,1)}(x)
      u_{k+1}(x) = Lu_k(x), k=0,1,2,...
      
      with the operator
      
      Lf(x) = 2*int_{-\infty}^x (f(2*t+1)-f(2*t-1)) dt
            = int_{-\infty}^{2*x} (f(t+1)-f(t-1)) dt
	    
    converges to a fixed point of L which is a C^\infty function up(x) on [-1,1].
    By definition, u_k is a spline of order k+1 with nodes at the dyadic points
    -1, -1+2^{-k+1}, -1+2*2^{-k+1},..., 1
 
    The following class provides a function object for u_k.
    
    References:
      Rvachev, "Compactly supported solutions of functional-differential equations",
      Russ. Math. Surv. 45 No. 1 (1990), 87-120
  */
  class ApproximateUpFunction
    : public Function<1,double>
  {
  public:
    ApproximateUpFunction(const unsigned int k);
    virtual ~ApproximateUpFunction();
    double value(const Point<1,double>& p,
		 const unsigned int component = 0) const;
    void vector_value(const Point<1,double> &p,
		      Vector<double>& values) const;
    
  private:
    unsigned int k_;
    Array1D<double> coeffs;
  };

}  

// include implementation
#include <numerics/up_function.cpp>
  
#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_BIHARMONIC_1D_TESTCASE_H
#define _FRAMETL_BIHARMONIC_1D_TESTCASE_H

#include <utils/function.h>
#include <functional.h>
#include <aggregated_frame.h>
#include <geometry/point.h>

using MathTL::Function;
using MathTL::Point;


namespace FrameTL
{
  /*!
    exact solution for the test case for the one-dimensional biharmonic equation
    u(x) = -cos(2*Pi*x)+1 + p(x)
    with
    p(x) = 48*x^4-63*x^3+47/2*x^2, for 0 <= x < 1/2,
    p(x) = 48*x^4-127*x^3+235/2*x^2-46*x+15/2, for 1/2 <= x <= 1
   */
  class Biharmonic1D_Solution
    : public Function<1>
  {
  public:
    //! default constructor
    Biharmonic1D_Solution()
    {
    }

    //! point evaluation
    double value(const Point<1>& p, const unsigned int component = 0) const;

    //! dummy concretisation for the compiler
    void vector_value(const Point<1> &p, Vector<double>& values) const
    { ; }
  };


  /*! A class inheriting from Function<1>, modeling the function f of the
    integrand in \int v(x) f(x) dx
   */
  class Biharmonic1D_RHS_Integrand
    : public Function<1>
  {
  public:
    //! point evaluation
    double value(const Point<1>& p, const unsigned int component = 0) const;

    //! dummy concretisation for the compiler
    void vector_value(const Point<1> &p, Vector<double>& values) const
    { ; }
  };


  /*!
    functional for the righthand side of the 1D test case for the biharmonic equation
    v -> 4*v'(1/2) - 384*v(1/2) - 16 * \int_0^1 (cos(2*Pi*x)*Pi^4-72) * v(x) dx
   */
  template<class IBASIS>
  class Biharmonic1D_RHS
    : public Functional<IBASIS, 1, 1>
  {
  public:
    /*!
      Constructor with an instance of a frame
     */
    Biharmonic1D_RHS(const AggregatedFrame<IBASIS, 1>* frame)
      : Functional<IBASIS,1,1>(new Biharmonic1D_RHS_Integrand(), frame)
    {}

    /*!
      Evaluate the functional with the wavelet with given index lambda
     */
    virtual double evaluate(const typename AggregatedFrame<IBASIS, 1>::Index& lambda) const;
  };
}

#include <biharmonic_1d_testcase.cpp> // include implementation

#endif // _FRAMETL_BIHARMONIC_1D_TESTCASE_H

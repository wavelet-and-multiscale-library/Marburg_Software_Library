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
    Exact solution for the test case for the one-dimensional biharmonic equation,
    \f$u(x) = -\cos(2\pi x)+1 + p(x)\f$
    with
    \f$p(x) = 48x^4-63x^3+47/2x^2\f$, for \f$0 \leq x < 1/2\f$,
    \f$p(x) = 48x^4-127x^3+235/2x^2-46x+15/2\f$, for \f$1/2 \leq x \leq 1\f$
   */
  class Biharmonic1D_Solution
    : public Function<1>
  {
  public:
    //! Default constructor.
    Biharmonic1D_Solution()
    {
    }

    //! Point evaluation.
    double value(const Point<1>& p, const unsigned int component = 0) const;

    //! Dummy concretisation for the compiler.
    void vector_value(const Point<1> &p, Vector<double>& values) const
    { ; }
  };


  /*! A class inheriting from Function<1>, modeling the function f of the
    integrand in \f$\int v(x) f(x) dx\f$.
   */
  class Biharmonic1D_RHS_Integrand
    : public Function<1>
  {
  public:
    //! Point evaluation.
    double value(const Point<1>& p, const unsigned int component = 0) const;

    //! Dummy concretisation for the compiler.
    void vector_value(const Point<1> &p, Vector<double>& values) const
    { ; }
  };


  /*!
    Functional for the righthand side of the 1D test case for the biharmonic equation
    \f$ v \mapsto 4v'(1/2) - 384v(1/2) - 16  \int_0^1 (\cos(2\pi x)\pi^4-72) * v(x) dx\f$.
   */
  template<class IBASIS>
  class Biharmonic1D_RHS
    : public Functional<IBASIS, 1, 1>
  {
  public:
    /*!
      Constructor with an instance of a frame.
     */
    Biharmonic1D_RHS(const AggregatedFrame<IBASIS, 1>* frame)
      : Functional<IBASIS,1,1>(new Biharmonic1D_RHS_Integrand(), frame)
    {}

    /*!
      Virtual destructor.
     */
    virtual ~Biharmonic1D_RHS()
    {}

    /*!
      Evaluate the functional with the wavelet with given index lambda.
     */
    virtual double evaluate(const typename AggregatedFrame<IBASIS, 1>::Index& lambda) const;
  };
}

#include <biharmonic_1d_testcase.cpp> // include implementation

#endif // _FRAMETL_BIHARMONIC_1D_TESTCASE_H

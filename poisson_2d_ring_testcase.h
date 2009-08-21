// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of WaveletTL - the Wavelet Template Library      |
// |                                                                    |
// | Copyright (c) 2002-2007                                            |
// | Thorsten Raasch, Manuel Werner, Andreas Schneider                  |
// +--------------------------------------------------------------------+

#ifndef _FRAMETL_POISSON_2D_RING_TESTCASE_H
#define _FRAMETL_POISSON_2D_RING_TESTCASE_H

#include <utils/function.h>
#include <aggregated_frame.h>
#include <geometry/point.h>
#include <numerics/corner_singularity.h>

using MathTL::Function;
using MathTL::Point;


namespace FrameTL
{
  class Poisson_Solution_Ring
    : public Function<2>
  {
  public:
    //! default constructor
    Poisson_Solution_Ring()
    {
    }
    //! point evaluation
    double value(const Point<2>& p, const unsigned int component = 0) const;
    
    //! dummy concretisation for the compiler
    void vector_value(const Point<2> &p, Vector<double>& values) const
    { ; }
  };

  class Poisson_RHS_Ring
    : public Function<2>
  {
  public:
    //! default constructor
    Poisson_RHS_Ring()
    {
    }
    //! point evaluation
    double value(const Point<2>& p, const unsigned int component = 0) const;
    
    //! dummy concretisation for the compiler
    void vector_value(const Point<2> &p, Vector<double>& values) const
    { ; }
  };

  class Poisson_SolutionGradient_Ring
    : public Function<2>
  {
  public:
    //! default constructor
    Poisson_SolutionGradient_Ring()
    {
    }
    void vector_value(const Point<2> &p, Vector<double>& values) const;

    //! point evaluation
    double value(const Point<2>& p, const unsigned int component = 0) const
    {
      return 0.0;
    }
    
    
  };



}
#include <poisson_2d_ring_testcase.cpp> // include implementation

#endif

// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of FrameTL - the Wavelet Template Library        |
// |                                                                    |
// | Copyright (c) 2002-2009                                            |
// | Thorsten Raasch, Manuel Werner                                     |
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

  /*!
    This class models the exact solution of a test problem for the Poisson
    equation in the two-dimensional ring-shaped domain
    \f$(-1,2)^2\setminus [0,1]^2\f$; cf. Sect. 7.2.3 in Manuel's PhD thesis.
    The solution consists of four singularity functions
    (as in eq. (4.4.4) in Manuel's PhD thesis) situated around the re-entrant corners.
    The whole funtion is plotted in Fig.7.22 (left).
  */
  class Poisson_Solution_Ring
    : public Function<2>
  {
  public:
    
    //! Default constructor.
    Poisson_Solution_Ring()
    {
    }

    //! Point evaluation.
    double value(const Point<2>& p, const unsigned int component = 0) const;
    
    //! Dummy concretization for the compiler.
    void vector_value(const Point<2> &p, Vector<double>& values) const
    { ; }
  };

  /*!
    This class models the right-hand side of a test problem for the Poisson
    equation in the two-dimensional ring-shaped domain
    \f$(-1,2)^2\setminus [0,1]^2\f$; cf. Sect. 7.2.3 in Manuel's PhD thesis.
    The whole funtion is plotted in Fig.7.22 (right).
  */
  class Poisson_RHS_Ring
    : public Function<2>
  {
  public:


    //! Default constructor.
    Poisson_RHS_Ring()
    {
    }

    //! Point evaluation.
    double value(const Point<2>& p, const unsigned int component = 0) const;
    
    //! Dummy concretization for the compiler.
    void vector_value(const Point<2> &p, Vector<double>& values) const
    { ; }
  };

  /*!
    This class models the gradient of the exact solution of a test
    problem for the Poisson equation in the two-dimensional ring-shaped domain
    \f$(-1,2)^2\setminus [0,1]^2\f$; cf. Sect. 7.2.3 in Manuel's PhD thesis.
  */
  class Poisson_SolutionGradient_Ring
    : public Function<2>
  {
  public:
    
    //! Default constructor.
    Poisson_SolutionGradient_Ring()
    {
    }

    //! Point evaluation.
    void vector_value(const Point<2>& p, Vector<double>& values) const;

    //! Dummy concretization for the compiler.
    double value(const Point<2>& p, const unsigned int component = 0) const
    {
      return 0.0;
    }
  };
}
#include <poisson_2d_ring_testcase.cpp> // include implementation

#endif

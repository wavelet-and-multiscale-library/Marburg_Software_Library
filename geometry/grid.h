// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2004                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_GRID_H
#define _MATHTL_GRID_H

#include <utils/array1d.h>
#include <geometry/point.h>

namespace MathTL
{
  /*!
    Abstract base class for n-dimensional rectangular grids in the
    style of Matlab. There, rectangular grids are the core ingredients for plotting
    functions from n-D to m-D.
    A 1-dimensional grid is just a vector x holding the mesh points.
    A 2-dimensional grid (a so-called quad-mesh) consists of 2 matrices x and y,
    holding the x- and y-coordinates of the mesh points.
    
    reference: Matlab/Octave help for the command 'surf'
   */
  template <unsigned int DIM>
  class Grid
  {
  public:
    /*!
      default constructor: empty grid
    */
    Grid();

    /*!
      Matlab output of the grid onto a stream
      
    */
    void matlab_output(std::ostream& os) const;
  };

  /*!
    specialization of Grid to one space dimension:
    1-dimensional grids are just vectors holding the mesh points.
  */
  template <>
  class Grid<1>
  {
  public:
    /*!
      default constructor: empty grid
    */
    Grid();

    /*!
      construct a grid from an array of 1D points
    */
    Grid(const Array1D<Point<1> >& grid);

    /*!
      Matlab output of the grid onto a stream
    */
    void matlab_output(std::ostream& os) const;
    
  private:
    /*!
      internal storage for the grid points
    */
    Array1D<Point<1> > grid_;
  };
}

// include implementation of inline functions
#include <geometry/grid.cpp>

#endif

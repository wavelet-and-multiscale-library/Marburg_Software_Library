// -*- c++ -*-

// +--------------------------------------------------------------------+
// | This file is part of MathTL - the Mathematical Template Library    |
// |                                                                    |
// | Copyright (c) 2002-2005                                            |
// | Thorsten Raasch                                                    |
// +--------------------------------------------------------------------+

#ifndef _MATHTL_GRID_H
#define _MATHTL_GRID_H

#include <geometry/point.h>
#include <utils/array1d.h>
#include <algebra/matrix.h>
#include <io/matrix_io.h>

namespace MathTL
{
  /*!
    Base class for n-dimensional rectangular grids in the
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
    Grid() {}

    /*!
      number of grid points
    */
    unsigned int size() const;

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
    Grid() : grid_() {}

    /*!
      construct a 1D grid from an array of 1D points
    */
    Grid(const Array1D<double>& grid) : grid_(grid) {}

    /*!
      construct an equidistant 1D grid with N+1 points
    */
    Grid(const double a, const double b, const unsigned N)
      : grid_(N+1)
    {
      for (unsigned int n(0); n <= N; n++)
	grid_[n] = a+(b-a)*n/N;
    }

    /*!
      number of grid points
    */
    inline unsigned int size() const { return grid_.size(); }

    /*!
      Matlab output of the grid onto a stream
    */
    void matlab_output(std::ostream& os) const
    {
      os << "x = "
	 << grid_
	 << ";"
	 << std::endl;
    }
    
  protected:
    /*!
      internal storage for the grid points
    */
    Array1D<double> grid_;
  };

  /*!
    specialization of Grid to two space dimensions:
    2-dimensional grids (quad-meshes) consist of 2 matrices x and y,
    holding the x- and y-coordinates of the mesh points
  */
  template <>
  class Grid<2>
  {
  public:
    /*!
      default constructor: empty grid
    */
    Grid() : gridx_(), gridy_() {}

    /*!
      construct a 2D grid from two matrices 
    */
    Grid(const Matrix<double>& gridx, const Matrix<double>& gridy)
      : gridx_(gridx), gridy_(gridy)
    {
    }

    /*!
      construct an equidistant 2D grid with (N_x+1)*(N_y+1) points
    */
    Grid(const Point<2>& a, const Point<2>& b,
	 const unsigned N_x, const unsigned N_y)
      : gridx_(N_y+1, N_x+1), gridy_(N_y+1, N_x+1)
    {
      for (unsigned int n_x(0); n_x <= N_x; n_x++)
	for (unsigned int n_y(0); n_y <= N_y; n_y++)
	  {
	    gridx_(n_y, n_x) = a(0) + (b(0)-a(0))*n_x/N_x;
	    gridy_(n_y, n_x) = a(1) + (b(1)-a(1))*n_y/N_y;
	  }
    }

    /*!
      number of grid points
    */
    inline unsigned int size() const { return gridx_.size() * gridy_.size(); }

    /*!
      Matlab output of the grid onto a stream
    */
    void matlab_output(std::ostream& os) const
    {
      os << "x = ";
      print_matrix(gridx_, os);
      os << ";" << std::endl;
      
      os << "y = ";
      print_matrix(gridy_, os);
      os << ";" << std::endl;
    }
    
  protected:
    /*!
      internal storage for the grid points
      (a 2D array would be sufficient, which is not available at the moment)
    */
    Matrix<double> gridx_, gridy_;
  };
}

#endif

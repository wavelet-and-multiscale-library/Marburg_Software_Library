// implementation of grid.h inline functions

#include <io/matrix_io.h>

namespace MathTL
{
  Grid<1>::Grid()
    : grid_()
  {
  }

  Grid<1>::Grid(const Array1D<double>& grid)
    : grid_(grid)
  {
  }

  Grid<1>::Grid(const double a, const double b, const unsigned int N)
    : grid_(N+1)
  {
    for (unsigned int n(0); n <= N; n++)
      grid_[n] = a+(b-a)*n/N;
  }

  Grid<2>::Grid()
    : gridx_(), gridy_()
  {
  }

  Grid<2>::Grid(const Matrix<double>& gridx,
		const Matrix<double>& gridy)
    : gridx_(gridx), gridy_(gridy)
  {
  }

  Grid<2>::Grid(const Point<2>& a, const Point<2>& b,
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

  inline unsigned int Grid<1>::size() const
  {
    return grid_.size();
  }

  inline unsigned int Grid<2>::size() const
  {
    return gridx_.size();
  }

  void Grid<1>::matlab_output(std::ostream& os) const
  {
    os << "x = "
       << grid_
       << ";"
       << std::endl;
  }

  void Grid<2>::matlab_output(std::ostream& os) const
  {
    os << "x = ";
    print_matrix(gridx_, os);
    os << ";" << std::endl;

    os << "y = ";
    print_matrix(gridy_, os);
    os << ";" << std::endl;
  }
}

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

  Grid<2>::Grid()
    : gridx_(), gridy_()
  {
  }

  Grid<2>::Grid(const Matrix<double>& gridx,
		const Matrix<double>& gridy)
    : gridx_(gridx), gridy_(gridy)
  {
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

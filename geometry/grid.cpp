// implementation of grid.h inline functions

namespace MathTL
{
  Grid<1>::Grid()
    : grid_()
  {
  }

  Grid<1>::Grid(const Array1D<Point<1> >& grid)
    : grid_(grid)
  {
  }

  void Grid<1>::matlab_output(std::ostream& os) const
  {
    os << "x = "
       << grid_
       << ";"
       << std::endl;
  }
}

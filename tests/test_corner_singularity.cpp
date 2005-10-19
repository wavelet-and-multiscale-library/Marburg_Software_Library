#include <iostream>
#include <fstream>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>
#include <numerics/corner_singularity.h>

using namespace MathTL;
using namespace std;

int main()
{
  cout << "Testing CornerSingularity ..." << endl;

  Point<2> origin;
  CornerSingularity s(origin, 0.5, 1.5);

  cout << "- Matlab output of the corner singularity on a 2D grid..." << endl;
  Grid<2> grid(Point<2>(-1.0, -1.0), Point<2>(1.0, 1.0), 1<<6, 1<<6);
  SampledMapping<2> h(grid, s);
  std::ofstream fs("corner.m");
  h.matlab_output(fs);
  fs << "surf(x,y,z)" << endl;
  fs.close();
  cout << "  ...done, see file corner.m!" << endl;

  return 0;
}

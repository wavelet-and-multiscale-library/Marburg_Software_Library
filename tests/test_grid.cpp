#include <iostream>
#include <utils/array1d.h>
#include <geometry/point.h>
#include <geometry/grid.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::Grid ..." << endl;

  cout << "- empty 1D grid:" << endl;
  Grid<1>().matlab_output(cout);
  
  cout << "- a non-empty 1D grid:" << endl;
  Array1D<Point<1> > points(3);
  points[0] = Point<1>(1.0);
  points[1] = Point<1>(2.4);
  points[2] = Point<1>(3.0);
  Grid<1>(points).matlab_output(cout);

  return 0;
}

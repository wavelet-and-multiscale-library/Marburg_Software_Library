#include <iostream>
#include <geometry/point.h>
#include <utils/array1d.h>
#include <utils/function.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::Function ..." << endl;
  cout << "- ZeroFunction<1> at x=23.0: "
       << ZeroFunction<1>().value(Point<1>(23.0)) << endl;
  cout << "- ZeroFunction<2> at (x,y)=(42.0, -1.5): "
       << ZeroFunction<2>().value(Point<2>(42.0, -1.5)) << endl;
  Array1D<double> y(2);
  ZeroFunction<2>(2).vector_value(Point<2>(42.0, -1.5), y);
  cout << "- vector-valued ZeroFunction<2> at (x,y)=(42.0, -1.5): "
       << y << endl;
}

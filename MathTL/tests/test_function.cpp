#include <iostream>
#include <geometry/point.h>
#include <algebra/vector.h>
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
  Vector<double> y(2);
  ZeroFunction<2>(2).vector_value(Point<2>(42.0, -1.5), y);
  cout << "- vector-valued ZeroFunction<2> at (x,y)=(42.0, -1.5): "
       << y << endl;

  Vector<double> z(2, "1 3.5");
  ConstantFunction<2> cf(z);
  cf.vector_value(Point<2>(42.0, -1.5), y);
  cout << "- vector-valued constant function at (x,y)=(42.0, -1.5): "
       << y << endl;
  
  
  Vector<double> z0(1, "3");
  ConstantFunction<2> cf0(z0);
  ProductFunction<2> pf0(&cf0, &cf0);
  cout << "- a product of two constants, evaluated at some point: "
       << pf0.value(Point<2>(1,2)) << endl;
}

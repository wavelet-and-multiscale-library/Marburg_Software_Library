#include <iostream>
#include <geometry/point.h>
#include <utils/function.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::Function ..." << endl;
  cout << "- ZeroFunction<1> at x=23.0: "
       << ZeroFunction<1>().value(Point<1>(23.0)) << endl;
  cout << "- ZeroFunction<2> st (x,y)=(42.0, -1.5): "
       << ZeroFunction<2>().value(Point<2>(42.0, -1.5)) << endl;
}

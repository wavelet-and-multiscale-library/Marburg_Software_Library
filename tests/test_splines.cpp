#include <iostream>
#include <geometry/point.h>
#include <numerics/splines.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::Spline ..." << endl;

  cout << "* some point values of the default constant spline:" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  Spline<1>().value(Point<1>(x)) << endl;
  }

  cout << "* some point values of the default linear spline:" << endl;
  for (double x = -0.5; x <= 2.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  Spline<2>().value(Point<1>(x)) << endl;
  }

  cout << "* some point values of the default quadratic spline:" << endl;
  for (double x = -0.5; x <= 3.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  Spline<3>().value(Point<1>(x)) << endl;
  }
}

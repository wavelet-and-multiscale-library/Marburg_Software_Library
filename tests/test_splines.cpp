#include <iostream>
#include <fstream>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>
#include <numerics/splines.h>
#include <utils/array1d.h>

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

  cout << "* a linear spline from a given knot sequence:" << endl;
  Array1D<double> knots(4+2);
  knots[0] =  0.0;
  knots[1] =  1.0;
  knots[2] =  1.5;
  knots[3] =  2.0;
  knots[4] =  2.5;
  knots[5] =  3.0;
  Array1D<double> coeffs(4);
  coeffs[0] = 1;
  coeffs[1] = 0;
  coeffs[2] = 0;
  coeffs[3] = 0;
  SampledMapping<1>(Grid<1>(0.0, 3.0, 30), Spline<2>(knots, coeffs)).matlab_output(cout);

  cout << "* writing point values of a linear spline with multiple knots to a file..." << endl;
  knots[0] =  0.0;
  knots[1] =  0.0;
  knots[2] =  1.0;
  knots[3] =  2.0;
  knots[4] =  2.0;
  knots[5] =  3.0;
  coeffs[0] = 0;
  coeffs[1] = 0;
  coeffs[2] = 0;
  coeffs[3] = 1;
  std::ofstream fs("multipleknots2.m");
  SampledMapping<1>(Grid<1>(0.0, 3.0, 60), Spline<2>(knots, coeffs)).matlab_output(fs);
  fs.close();

  cout << "* writing point values of a cubic spline with multiple knots to a file..." << endl;
  knots.resize(7);
  knots[0] =  0.0;
  knots[1] =  0.0;
  knots[2] =  0.0;
  knots[3] =  1.0;
  knots[4] =  2.0;
  knots[5] =  2.0;
  knots[6] =  2.0;
  coeffs.resize(4);
  coeffs[0] = 1;
  coeffs[1] = 1;
  coeffs[2] = 2;
  coeffs[3] = 1;
  std::ofstream fs2("multipleknots3.m");
  SampledMapping<1>(Grid<1>(0.0, 2.0, 1000), Spline<3>(knots, coeffs)).matlab_output(fs2);
  fs2.close();
}

#include <iostream>
#include <geometry/chart.h>
#include <geometry/point.h>
#include <algebra/matrix.h>
#include <algebra/vector.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing some charts:" << endl;

  AffineLinearMapping<2> kappa_0;
  cout << "* the identity mapping: " << endl << kappa_0;
  Point<2> x(1.0, 2.0), y;
  kappa_0.map_point(x, y);
  cout << "  x=" << x << " is mapped to y=" << y << endl;
  kappa_0.map_point_inv(y, x);
  cout << "  y=" << y << " is mapped back to x=" << x << endl;

  Matrix<double> A(2, 2, "2 0 1 -4");
  Point<2> b(0, 1);
  AffineLinearMapping<2> kappa_1(A, b);
  cout << "* affine linear map: " << endl << kappa_1;
  kappa_1.map_point(x, y);
  cout << "  x=" << x << " is mapped to y=" << y << endl;
  kappa_1.map_point_inv(y, x);
  cout << "  y=" << y << " is mapped back to x=" << x << endl;
  
  Point<2> P(1.0, 2.0);
  kappa_1.map_point_inv(P, y);
  cout << "  test point: P=" << P
       << " (inverse: " << y << ")"
       << ", in patch: " << kappa_1.in_patch(P) << endl;

  // [0,1] -> [-1,1]
  Matrix<double> C(1, 1, "2");
  Point<1> d(-1.0);
  AffineLinearMapping<1> kappa_2(C, d);
  cout << "* affine linear map with C=" << endl << C
       << "  and d=" << d << ":" << endl;
  Point<1> x2(0.5), y2;
  kappa_2.map_point(x2, y2);
  cout << "  x=" << x2 << " is mapped to y=" << y2 << endl;
  kappa_2.map_point_inv(y2, x2);
  cout << "  y=" << y2 << " is mapped back to x=" << x2 << endl;

  return 0;
}

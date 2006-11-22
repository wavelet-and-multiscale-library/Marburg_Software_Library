#include <iostream>
#include <algebra/matrix.h>
#include <algebra/vector.h>
#include <algebra/qs_matrix.h>

using namespace std;
using namespace MathTL;

int main()
{
  cout << "Constructing a quasi-stationary matrix..." << endl;
  
  Matrix<double> one(1, 1, "1");
  Matrix<double> two(1, 1, "2");
  Vector<double> b1(2, "3 4");
  Vector<double> b2(2, "5 6");

  QuasiStationaryMatrix<double> Q(2, 13, 5, one, two, b1, b2, 2, 1);
  cout << "Q=" << endl << Q;

  Q.set_level(3);
  cout << "Q on next higher level:" << endl << Q;

  return 0;
}

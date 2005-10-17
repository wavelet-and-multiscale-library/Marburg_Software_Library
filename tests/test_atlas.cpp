#include <iostream>
#include <geometry/chart.h>
#include <geometry/atlas.h>
#include <geometry/point.h>
#include <algebra/matrix.h>
#include <algebra/symmetric_matrix.h>
#include <utils/array1d.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing the atlas class:" << endl;

  Atlas<2> E;
  cout << "* the empty atlas:" << endl << E;

  AffineLinearMapping<2> I;
  Array1D<Chart<2>*> kappas(1);
  kappas[0] = &I;
  SymmetricMatrix<bool> S(1);
  S(0, 0) = true;
  Atlas<2> A(kappas, S);
  cout << "* the identity atlas:" << endl << A;

  return 0;
}

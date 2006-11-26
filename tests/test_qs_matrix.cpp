#include <iostream>
#include <algebra/matrix.h>
#include <algebra/vector.h>
#include <algebra/qs_matrix.h>
#include <algebra/sparse_matrix.h>

using namespace std;
using namespace MathTL;

int main()
{
  cout << "Constructing a quasi-stationary matrix..." << endl;
  
  Matrix<double> one(3, 2, "1 2 3 4 5 6");
  Matrix<double> two(4, 1, "-4 -3 -2 -1");
  Vector<double> b1(3, "3 4 5");
  Vector<double> b2(3, "5 6 -7");

  QuasiStationaryMatrix<double> Q(3, 17, 9, one, two, b1, b2, 2, 1);
  cout << "Q=" << endl << Q;

//   Q.set_level(4);
//   cout << "Q on next higher level:" << endl << Q;

  SparseMatrix<double> S;
  Q.to_sparse(S);
  cout << "Q as a sparse matrix:" << endl << S;

  Vector<double> x(Q.column_dimension());
  x[4] = x[5] = 1;
  Vector<double> Qx(Q.row_dimension());
  Q.apply(x, Qx);
  cout << "Q applied to x=" << x << " yields Qx=" << Qx << endl;

  Vector<double> y(Q.row_dimension());
  y[11] = y[12] = 1;
  Vector<double> Qty(Q.column_dimension());
  Q.apply_transposed(y, Qty);
  cout << "Q^T applied to y=" << y << " yields Q^Ty=" << Qty << endl;

  x.resize(Q.column_dimension()+2);
  x[5] = x[6] = 1;
  Qx.resize(Q.row_dimension()+2);
  Q.apply(x, Qx, 1, 2);
  cout << "Q applied to the part of x=" << x << " from index 1," << endl
       << "written into Qx from index 2 on yields Qx=" << Qx << endl;

  y.resize(Q.row_dimension()+3);
  y[14] = y[15] = 1;
  Qty.resize(Q.column_dimension()+1);
  Q.apply_transposed(y, Qty, 3, 1);
  cout << "Q^T applied to the part of y=" << y << " from index 3," << endl
       << "written into Qty from index 1 on yields Qty=" << Qty << endl;

  return 0;
}

#include <iostream>
#include <numerics/decomposable_matrix.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing decomposable matrices ..." << endl;

  DecomposableMatrix<double> N(2);
  cout << "- a default decomposable matrix N:" << endl
       << N;

  DecomposableMatrix<double> A(Matrix<double>(4, 3, "0 1 2 3 4 5 6 7 9 0 1 0"));
  cout << "- a rectangular matrix A:" << endl
       << A;

  A.decompose(DecomposableMatrix<double>::LU);
  A.decompose(DecomposableMatrix<double>::none);
  
  cout << "- after LU decomposition and reversion:" << endl
       << A;

  DecomposableMatrix<double> M(Matrix<double>(3, 3, "2 -1 0 -1 2 -1 0 -1 2"));
  cout << "- the standard matrix M_3:" << endl
       << M;
  Vector<double> b(3, "3 0 -1"), x;
  M.solve(b, x);
  cout << "- solving Mx=" << b << " yields x=" << x << endl; 

  DecomposableMatrix<double> B(Matrix<double>(3, 3, "0 -1 2 2 -1 0 -1 2 -1"));
  cout << "- the standard matrix B=PM_3:" << endl
       << B;
  b = Vector<double>(3, "-1 3 0");
  B.solve(b, x);
  cout << "- solving Bx=" << b << " yields x=" << x << endl; 

  A = Matrix<double>(4, 3, "2 2 1 0 0 1 1 -1 0 -2 0 1");
  cout << "- a rectangular matrix A:" << endl
       << A;

  A.decompose(DecomposableMatrix<double>::QU);
  A.decompose(DecomposableMatrix<double>::none);

  cout << "- after QU decomposition and reversion:" << endl
       << A;

  cout << "- the standard matrix again:" << endl
       << M;

  M.decompose(DecomposableMatrix<double>::QU);
  M.decompose(DecomposableMatrix<double>::none);

  cout << "- after QU decomposition and reversion:" << endl
       << M;

  b = Vector<double>(3, "3 0 -1");
  M.solve(b, x);
  cout << "- solving Mx=" << b << " with QU decomposition yields x=" << x << endl; 

  return 0;
}

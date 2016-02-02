#include <iostream>
#include <algebra/vector.h>
#include <algebra/tridiagonal_matrix.h>
#include <numerics/preconditioner.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing some preconditioners..." << endl;

  TridiagonalMatrix<double> S(4);
  S.set_entry(0, 0, 2.0);
  S.set_entry(1, 1, 2.0);
  S.set_entry(2, 2, 2.0);
  S.set_entry(3, 3, 2.0);
  S.set_entry(0, 1, -1.0);
  S.set_entry(1, 2, -1.0);
  S.set_entry(2, 3, -1.0);
  S.set_entry(1, 0, -1.0);
  S.set_entry(2, 1, -1.0);
  S.set_entry(3, 2, -1.0);
  cout << "- the standard matrix:" << endl << S;

  JacobiPreconditioner<TridiagonalMatrix<double>,Vector<double> > JP(S);
  Vector<double> x(4, "1 2 3 4"), y;
  JP.apply_preconditioner(x, y);
  cout << "- applying the Jacobi preconditioner to x=" << x << " yields y=P^{-1}x=" << y << endl;

  return 0;
}

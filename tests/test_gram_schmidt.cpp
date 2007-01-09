#include <iostream>
#include "numerics/gram_schmidt.h"
#include "algebra/vector.h"
#include "algebra/vector_norms.h"
#include "algebra/piecewise.h"

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the Gram-Schmidt process ..." << endl;

  cout << endl << "Applying the process on some vectors ..." << endl;
  Array1D< Vector<double> > v(2);
  v[0] = Vector<double>(2, "3 1");
  v[1] = Vector<double>(2, "2 2");

  cout << "v[0] = " << v[0] << ", v[1] = " << v[1] << endl;

  cout << "After orthogonalization ..." << endl;
  gramSchmidtProcess(v);
  // exact result: [3,1]/sqrt(10), [-1,3]/sqrt(10)
  cout << "v[0] = " << v[0] << ", v[1] = " << v[1] << endl;

  cout << endl << "Applying the process on some functions (polynomials on [-1,1]) ..." << endl;

  // monomial basis
  Polynomial<double> x0;
  x0.set_coefficient(0,1);
  Polynomial<double> x1;
  x1.set_coefficient(1,1);
  Polynomial<double> x2;
  x2.set_coefficient(2,1);

  // ... on [-1,1]
  Array1D< Piecewise<double> > b(3);
  b[0].set_local_expansion(-1, x0);
  b[0].set_local_expansion(0, x0);
  b[1].set_local_expansion(-1, x1);
  b[1].set_local_expansion(0, x1);
  b[2].set_local_expansion(-1, x2);
  b[2].set_local_expansion(0, x2);
  cout << b << endl;

  gramSchmidtProcess(b);
  cout << "After orthogonalization ..." << endl;
  cout << b << endl;

  return 0;
}

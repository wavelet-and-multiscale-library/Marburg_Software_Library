#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <algebra/fixed_matrix.h>
#include <algebra/fixed_vector.h>
#include <algebra/matrix_norms.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the fixed matrix class ..." << endl;

  cout << "- the default matrix M (empty):" << endl;
  FixedMatrix<double, 2> M;
  cout << M;

  cout << "- a 2x3 matrix N:" << endl;
  FixedMatrix<double, 2,3> N;
  cout << N;
  
  cout << "- writing single entries:" << endl;
  N(0,0) = 23; N(0,1) = 42; N(0,2) = 3; N(1,1) = -1.5;
  cout << N;

  cout << "- testing copy constructor:" << endl;
  FixedMatrix<double,2,3> Ncopy(N);
  cout << Ncopy;

  cout << "- are the two matrices equal?" << endl;
  if (N == Ncopy)
    cout << "  ... yes!" << endl;
  else
    cout << "  ... no!" << endl;

  FixedMatrix<double,2,3> O;
  O(0,0) = 3.14; O(1,1) = -7.0;
  cout << "- another matrix O:" << endl << O;
  N = O;
  cout << "- N=O:" << endl << N;

  O.scale(-2.0);
  cout << "- O scaled by -2.0:" << endl << O;

  FixedVector<double,3> x;
  FixedVector<double,2> y;
  x(0) = 1; x(1) = -1; x(2) = -2;
  cout << "- testing apply(), using the vector x="
       << x << ";" << endl;
  cout << "  Nx=";
  N.apply(x,y);
  cout << y << endl;
  
  cout << "- choose another matrix A (constructed from a string):" << endl;
  FixedMatrix<double,3,3> A("1 2 3 4 5 6 7 8 9");
  cout << A;
  cout << "  Ax=";
  FixedVector<double, 3> z;
  A.apply(x,z);
  cout << z << endl;
  cout << "  A^Tx=";
  A.apply_transposed(x,z);
  cout << z << endl;

  cout << "- A*A:" << endl
       << A*A << endl;

  cout << "- A^T:" << endl
       << transpose(A) << endl;

  return 0;
}

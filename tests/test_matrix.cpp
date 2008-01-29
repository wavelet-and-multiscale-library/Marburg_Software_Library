#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <algebra/matrix.h>
#include <algebra/triangular_matrix.h>
#include <algebra/vector.h>
#include <algebra/matrix_norms.h>
#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/tridiagonal_matrix.h>
#include <algebra/infinite_vector.h>
#include <algebra/kronecker_matrix.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the matrix classes ..." << endl;

  cout << "- the default matrix M (empty):" << endl;
  Matrix<double> M;
  cout << M;

  cout << "- memory consumption of M: " << M.memory_consumption() << endl;

  cout << "- a 2x3 matrix N:" << endl;
  Matrix<double> N(2,3);
  cout << N;
  
  cout << "- memory consumption of N: " << N.memory_consumption() << endl;

  cout << "- writing single entries:" << endl;
  N(0,0) = 23; N(0,1) = 42; N(0,2) = 3; N(1,1) = -1.5;
  cout << N;

  cout << "- testing copy constructor:" << endl;
  Matrix<double> Ncopy(N);
  cout << Ncopy;

  cout << "- are the two matrices equal?" << endl;
  if (N == Ncopy)
    cout << "  ... yes!" << endl;
  else
    cout << "  ... no!" << endl;

  Matrix<double> O(2,3);
  O(0,0) = 3.14; O(1,1) = -7.0;
  cout << "- another matrix O:" << endl << O;
  swap(N,O);
  cout << "- after swapping, N=" << endl << N
       << "  ... and O=" << endl << O;
  N = O;
  cout << "- N=O:" << endl << N;

  O.scale(-2.0);
  cout << "- O scaled by -2.0:" << endl << O;

  Matrix<double> Q;
  O.reflect(Q);
  cout << "- O reflected:" << endl << Q;

  Q.diagonal(4, -12.34);
  cout << "- Q is now a diagonal matrix:" << endl << Q;

  Vector<double> x(3), y(2);
  x(0) = 1; x(1) = -1; x(2) = -2;
  cout << "- testing apply(), using the vector x="
       << x << ";" << endl;
  cout << "  Nx=";
  N.apply(x,y);
  cout << y << endl;
  
  cout << "- choose another matrix A (constructed from a string):" << endl;
  Matrix<double> A(3,3,"1 2 3 4 5 6 7 8 9");
  cout << A;
  y.resize(3);
  cout << "  Ax=";
  A.apply(x,y);
  cout << y << endl;
  cout << "  A^Tx=";
  A.apply_transposed(x,y);
  cout << y << endl;

  cout << "- testing matrix norms for A:" << endl;
  cout << "  max. row sum ||A||_infty: " << row_sum_norm(A) << endl;
  cout << "  max. column sum ||A||_1: " << column_sum_norm(A) << endl;
  cout << "  Frobenius norm ||A||_F: " << frobenius_norm(A) << endl;

  cout << "- A*A:" << endl
       << A*A << endl;

  cout << "- A^T:" << endl
       << transpose(A) << endl;

// //   cout << "- testing row iterators:" << endl;
// //   for (RawMatrix<double>::row_const_iterator it(N.rowbegin());
// //        it != N.rowend(); ++it)
// //     {
// //       cout << "  row " << it.index() << ": ";
// //       for (RawMatrix<double>::row_const_iterator::entry_const_iterator ite(it.begin());
// // 	   ite != it.end(); ++ite)
// // 	cout << ite.entry() << " ";
// //       cout << endl;
// //     }
  
// //   cout << "- testing column iterators:" << endl;
// //   for (RawMatrix<double>::col_const_iterator it(N.colbegin());
// //        it != N.colend(); ++it)
// //     {
// //       cout << "  column " << it.index() << ": ";
// //       for (RawMatrix<double>::col_const_iterator::entry_const_iterator ite(it.begin());
// // 	   ite != it.end(); ++ite)
// // 	cout << ite.entry() << " ";
// //       cout << endl;
// //     }

  cout << "- a lower triangular default matrix:" << endl;
  LowerTriangularMatrix<double> Ldefault;
  cout << Ldefault;
  
  cout << "- a lower triangular 3x3 matrix L:" << endl;
  LowerTriangularMatrix<double> L(3);
  cout << L;

  cout << "- constructing a lower triangular matrix from a string:" << endl;
  LowerTriangularMatrix<double> Lbyrow(4,4,"1 2 3 4 5 6 7 8 9 10");
  cout << Lbyrow;

  cout << "- constructing a lower triangular matrix from a string (columnwise):" << endl;
  LowerTriangularMatrix<double> Lbycol(4,4,"1 2 4 7 3 5 8 6 9 10", false);
  cout << Lbycol;

  cout << "- testing apply(), using the vector ";
  L = Lbyrow;
  x.resize(L.column_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(L.row_dimension());
  L.apply(x,y);
  cout << y << endl;

  cout << "- testing apply_transposed(), result: ";
  y.resize(L.column_dimension());
  L.apply_transposed(x,y);
  cout << y << endl;

  cout << "- a lower triangular matrix, again:" << endl;
  cout << Lbyrow;
  LowerTriangularMatrix<double> Linv;
  Lbyrow.inverse(Linv);
  cout << "- the inverse of that matrix:" << endl;
  cout << Linv;

  cout << "- an upper triangular default matrix:" << endl;
  UpperTriangularMatrix<double> Rdefault;
  cout << Rdefault;
  
  cout << "- an upper triangular 3x3 matrix R:" << endl;
  UpperTriangularMatrix<double> R(3);
  cout << R;
  
  cout << "- constructing an upper triangular matrix from a string:" << endl;
  UpperTriangularMatrix<double> Rbyrow(4,4,"1 2 3 4 5 6 7 8 9 10");
  cout << Rbyrow;

  cout << "- constructing an upper triangular matrix from a string (columnwise):" << endl;
  UpperTriangularMatrix<double> Rbycol(4,4,"1 2 5 3 6 8 4 7 9 10", false);
  cout << Rbycol;

  cout << "- an upper triangular matrix, again:" << endl;
  cout << Rbyrow;
  UpperTriangularMatrix<double> Rinv;
  Rbyrow.inverse(Rinv);
  cout << "- the inverse of that matrix:" << endl;
  cout << Rinv;
  
  cout << "- testing apply(), using the vector ";
  R = Rbyrow;
  x.resize(R.column_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(R.row_dimension());
  R.apply(x,y);
  cout << y << endl;

  cout << "- testing apply_transposed(), result: ";
  y.resize(R.column_dimension());
  R.apply_transposed(x,y);
  cout << y << endl;

  cout << "- testing non-quadratic lower triangular matrices:" << endl;
  LowerTriangularMatrix<double> Lnonq1(3,4,"1 2 3 4 5 6");
  cout << Lnonq1;
  cout << "  (this has " << Lnonq1.size() << " nonzero entries)" << endl;
  LowerTriangularMatrix<double> Lnonq2(4,3,"1 2 3 4 5 6 7 8 9");
  cout << Lnonq2;
  cout << "  (this has " << Lnonq2.size() << " nonzero entries)" << endl;

  cout << "- first matrix, testing apply() for x=" << x << ", result: ";
  y.resize(Lnonq1.row_dimension());
  Lnonq1.apply(x,y);
  cout << y << endl;
  cout << "- first matrix, testing apply_transposed() for x=";
  x.resize(Lnonq1.row_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(Lnonq1.column_dimension());
  Lnonq1.apply_transposed(x,y);
  cout << y << endl;
  cout << "- second matrix, testing apply() for x=";
  x.resize(Lnonq2.column_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(Lnonq2.row_dimension());
  Lnonq2.apply(x,y);
  cout << y << endl;
  cout << "- second matrix, testing apply_transposed() for x=";
  x.resize(Lnonq2.row_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(Lnonq2.column_dimension());
  Lnonq2.apply_transposed(x,y);
  cout << y << endl;

  cout << "- testing non-quadratic upper triangular matrices:" << endl;
  UpperTriangularMatrix<double> Rnonq1(3,4,"1 2 3 4 5 6 7 8 9");
  cout << Rnonq1;
  cout << "  (this has " << Rnonq1.size() << " nonzero entries)" << endl;
  UpperTriangularMatrix<double> Rnonq2(4,3,"1 2 3 4 5 6");
  cout << Rnonq2;
  cout << "  (this has " << Rnonq2.size() << " nonzero entries)" << endl;

  cout << "- first matrix, testing apply() for x=";
  x.resize(Rnonq1.column_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(Rnonq1.row_dimension());
  Rnonq1.apply(x,y);
  cout << y << endl;
  cout << "- first matrix, testing apply_transposed() for x=";
  x.resize(Rnonq1.row_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(Rnonq1.column_dimension());
  Rnonq1.apply_transposed(x,y);
  cout << y << endl;
  cout << "- second matrix, testing apply() for x=";
  x.resize(Rnonq2.column_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(Rnonq2.row_dimension());
  Rnonq2.apply(x,y);
  cout << y << endl;
  cout << "- second matrix, testing apply_transposed() for x=";
  x.resize(Rnonq2.row_dimension());
  for (unsigned int i(0); i < x.size(); i++) x(i)=i+1;
  cout << x << ", result: ";
  y.resize(Rnonq2.column_dimension());
  Rnonq2.apply_transposed(x,y);
  cout << y << endl;

  cout << "- a symmetric default matrix:" << endl;
  SymmetricMatrix<double> S;
  cout << S;
  
  cout << "- a symmetric 3x3 matrix T:" << endl;
  SymmetricMatrix<double> T(3);
  cout << T;

  cout << "- constructing a symmetric matrix from a string:" << endl;
  SymmetricMatrix<double> Sbyrow(4,"1 2 3 4 5 6 7 8 9 10");
  cout << Sbyrow;

  cout << "- constructing a symmetric matrix from a string (columnwise):" << endl;
  SymmetricMatrix<double>
    Sbycol(3,"1 2 4 3 5 6",false); // note the permutation!
  cout << Sbycol;

  
  cout << "- a sparse default matrix:" << endl;
  SparseMatrix<double> W;
  cout << W;
  
  cout << "- an empty 2x2 matrix X:" << endl;
  SparseMatrix<double> X(2);
  cout << X;

  cout << "- set some entries to nontrivial values:" << endl;
  X.set_entry(0, 0, -23.0);
  X.set_entry(1, 0, 42.0);
  cout << X;

  cout <<"- test copy constructor:" << endl;
  SparseMatrix<double> Ycopy(X);
  cout << Ycopy;

  cout << "- after X.resize(3,2):" << endl;
  X.resize(3, 2);
  cout << X;

  cout << "- refill X again with some values:" << endl;
  X.set_entry(0, 0, -3);
  X.set_entry(0, 1, 2);
  X.set_entry(1, 0, 1);
  X.set_entry(1, 1, -1.5);
  X.set_entry(2, 1, 3);
  cout << X;

  x.resize(2);
  x[0] = 2; x[1] = 3;
  cout << "- applying X to x=" << x << " yields" << endl;
  X.apply(x, y);
  cout << y << endl;

  cout << "- applying X^T to y=" << y << " yields" << endl;
  X.apply_transposed(y, x);
  cout << x << endl;

  cout << "- a 2x2 sparse diagonal matrix:" << endl;
  SparseMatrix<double> Y;
  Y.diagonal(2, 42.0);
  cout << Y;

  cout << "- putting different subblocks into this matrix:" << endl;
  Matrix<double> sub1(1, 1, "23");
  SymmetricMatrix<double> sub2(1, "7");
  LowerTriangularMatrix<double> sub3(1, 1, "-1");
  Y.set_block(0, 1, sub1);
  Y.set_block(1, 0, sub2);
  Y.set_block(1, 1, sub3);
  cout << Y;

  cout << "- putting a block into this matrix in reflect mode:" << endl;
  Matrix<double> sub4(2, 2, "1 2 3 4");
  Y.set_block(0, 0, sub4, true);
  cout << Y;

  cout << "- a 3x2 sparse matrix F1:" << endl;
  SparseMatrix<double> F1(3, 2);
  F1.set_entry(0, 0, 1.0);
  F1.set_entry(1, 0, 2.0);
  F1.set_entry(1, 1, 3.0);
  F1.set_entry(2, 1, 4.0);
  cout << F1;

  cout << "- a 2x2 sparse matrix F2:" << endl;
  SparseMatrix<double> F2(2);
  F2.set_entry(0, 0, 1.5);
  F2.set_entry(1, 0, 1.0);
  F2.set_entry(1, 1, -2.5);
  cout << F2;

  SparseMatrix<double> F3;
  F3 = F1 * F2;
  cout << "- matrix product F1*F2:" << endl
       << F3;

  SparseMatrix<double> small(2, 2);
  small.set_entry(0, 0, 1e-5);
  small.set_entry(0, 1, 2e-5);
  small.set_entry(1, 0, 3e-5);
  small.set_entry(1, 1, -1e-5);
  cout << "- a sparse matrix with small entries:" << endl
       << small;

  InfiniteVector<double, Vector<double>::size_type> v;
  small.get_row(0, v);
  cout << "- extracted a row from small as a sparse vector:" << endl
       << v;

  small.get_row(1, v, 10);
  cout << "- another row, now with offset 10:" << endl << v;

  std::list<size_t> indices;
  std::list<double> entries;
  indices.push_back(1);
  entries.push_back(-2.7);
  small.set_row(0, indices, entries);
  cout << "- the matrix small after setting the first row:" << endl
       << small;

  small.compress(1.5e-5);
  cout << "- small after compressing with 1.5e-5:" << endl
       << small;

  small.get_row(1, v);
  cout << "- extract last row again:" << endl << v;

  Matrix<double> extract;
  small.get_block(0, 1, 2, 1, extract);
  cout << "- extracted a sub-block from small:" << endl
       << extract;

  cout << "- the matrix F2 again:" << endl
       << F2
       << "  and its transpose:" << endl
       << transpose(F2);
  
  cout << "- a tridiagonal default matrix:" << endl;
  TridiagonalMatrix<double> T1;
  cout << T1;

  cout << "- an empty 4x4 tridiagonal matrix T2:" << endl;
  TridiagonalMatrix<double> T2(4);
  cout << T2;

  cout << "- set some entries to nontrivial values:" << endl;
  T2.set_entry(0, 0, 1.0);
  T2.set_entry(1, 1, 2.0);
  T2.set_entry(2, 2, 3.0);
  T2.set_entry(3, 3, 4.0);
  T2.set_entry(0, 1, -1.0);
  T2.set_entry(1, 2, -2.0);
  T2.set_entry(2, 3, -3.0);
  T2.set_entry(1, 0, -3.5);
  T2.set_entry(2, 1, -2.5);
  T2.set_entry(3, 2, -1.5);
  cout << T2;

  cout <<"- test copy constructor:" << endl;
  TridiagonalMatrix<double> T2copy(T2);
  cout << T2copy;

  x.resize(4); y.resize(4);
  x[0] = x[1] = x[2] = x[3] = 1;
  cout << "- applying T2 to x=" << x << " yields" << endl;
  T2.apply(x, y);
  cout << y << endl;

  cout << "- applying T2^T to x=" << x << " yields" << endl;
  T2.apply_transposed(x, y);
  cout << y << endl;

  typedef Matrix<double> MATRIX;
  MATRIX M1(2, 3, "1 2 3 4 5 6"), M2(2, 2, "1 2 0 3");
  cout << "- a matrix M1=" << endl << M1;
  cout << "- a matrix M2=" << endl << M2;
  KroneckerMatrix<double,MATRIX,MATRIX> K(M1,M2);
  cout << "- Kronecker product K of M1 and M2:" << endl << K;

  x.resize(6); x[3] = 1;
  K.apply(x, y);
  cout << "- K applied to x=" << x << " yields Kx=" << y << endl;
  x.resize(6); x[4] = 1;
  K.apply(x, y);
  cout << "- K applied to x=" << x << " yields Kx=" << y << endl;

  y.resize(6);
  x.resize(4); x[0] = 1;
  K.apply_transposed(x, y);
  cout << "- K^T applied to x=" << x << " yields (K^T)x=" << y << endl;
  x.resize(4); x[1] = 1;
  K.apply_transposed(x, y);
  cout << "- K^T applied to x=" << x << " yields (K^T)x=" << y << endl;
  
  Matrix<double> KM(K);
  cout << "- Kronecker product of M1 and M2 as a Matrix<double> again:"
       << endl << KM;

  Vector<double> w(3, "1 2 3");
  Matrix<double> wMatrix(w);
  cout << "- construct a matrix from a vector:" << endl
       << wMatrix;
  wMatrix.reshape(1);
  cout << "- reshape this matrix:" << endl
       << wMatrix;

  M = Matrix<double>(2, 2, "1 3 2 4");
  cout << "- another matrix M=" << endl << M;
  Vector<double> colM; M.col(colM);
  cout << "- col(M)=" << colM << endl;
  M.resize(45,67);
  M.decol(colM, 2);
  cout << "- decol(col(M))=" << endl << M;

#if 1
  F1.resize(5,2);
  F1.set_entry(0,0,1);
  F1.set_entry(0,1,2);
  F1.set_entry(1,1,3);
  F1.set_entry(2,0,4);
  F1.set_entry(4,1,5);

  cout << "- a sparse matrix F1 again:" << endl << F1;
  cout << "- write F1 to a file..." << endl;
  F1.matlab_output("F1", "F1", 1);
  cout << "  ... done" << endl;
  F1.resize(1,1);
  cout << "- resized F1:" << endl << F1;
  cout << "- read F1 from the file again..." << endl;
  F1.matlab_input("F1");
  cout << "  ...done, F1=" << endl << F1;
#endif

  return 0;
}

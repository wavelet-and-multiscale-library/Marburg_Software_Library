#include <iostream>
#include <algebra/matrix.h>
#include <algebra/vector.h>
// #include "symmarray2d.h"
// #include "triangarray2d.h"

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

//   cout << "- A*A:" << endl
//        << A*A << endl;

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

//   cout << "- a symmetric default matrix:" << endl;
//   RawMatrix<double,SymmetricArray2D<double> > S;
//   cout << S << endl;
  
//   cout << "- a symmetric 3x3 matrix T:" << endl;
//   RawMatrix<double,SymmetricArray2D<double> > T(3,3);
//   cout << T << endl;

//   cout << "- constructing a symmetric matrix from a string:" << endl;
//   RawMatrix<double,SymmetricArray2D<double> > Sbyrow(4,4,"1 2 3 4 5 6 7 8 9 10");
//   cout << Sbyrow << endl;

//   cout << "- constructing a symmetric matrix from a string (columnwise):" << endl;
//   RawMatrix<double,SymmetricArray2D<double> >
//     Sbycol(3,3,"1 2 4 3 5 6",false); // note the permutation!
//   cout << Sbycol << endl;
  
//   cout << "- a lower triangular default matrix:" << endl;
//   RawMatrix<double,LowerTriangularArray2D<double> > Ldefault;
//   cout << Ldefault << endl;
  
//   cout << "- a lower triangular 3x3 matrix L:" << endl;
//   RawMatrix<double,LowerTriangularArray2D<double> > L(3,3);
//   cout << L << endl;

//   cout << "- constructing a lower triangular matrix from a string:" << endl;
//   RawMatrix<double,LowerTriangularArray2D<double> > Lbyrow(4,4,"1 2 3 4 5 6 7 8 9 10");
//   cout << Lbyrow << endl;

//   cout << "- constructing a lower triangular matrix from a string (columnwise):" << endl;
//   RawMatrix<double,LowerTriangularArray2D<double> > Lbycol(4,4,"1 2 4 7 3 5 8 6 9 10",false);
//   cout << Lbycol << endl;

//   cout << "- testing APPLY, using the vector ";
//   L = Lbyrow;
//   x.resize(L.coldim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   L.APPLY(x,y);
//   cout << y << endl;

//   cout << "- testing APPLYtr, result: ";
//   L.APPLYtr(x,y);
//   cout << y << endl;

//   cout << "- an upper triangular default matrix:" << endl;
//   RawMatrix<double,UpperTriangularArray2D<double> > Rdefault;
//   cout << Rdefault << endl;
  
//   cout << "- an upper triangular 3x3 matrix R:" << endl;
//   RawMatrix<double,UpperTriangularArray2D<double> > R(3,3);
//   cout << R << endl;
  
//   cout << "- constructing an upper triangular matrix from a string:" << endl;
//   RawMatrix<double,UpperTriangularArray2D<double> > Rbyrow(4,4,"1 2 3 4 5 6 7 8 9 10");
//   cout << Rbyrow << endl;

//   cout << "- constructing an upper triangular matrix from a string (columnwise):" << endl;
//   RawMatrix<double,UpperTriangularArray2D<double> > Rbycol(4,4,"1 2 5 3 6 8 4 7 9 10",false);
//   cout << Rbycol << endl;

//   cout << "- testing APPLY, using the vector ";
//   R = Rbyrow;
//   x.resize(R.coldim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   R.APPLY(x,y);
//   cout << y << endl;

//   cout << "- testing APPLYtr, result: ";
//   R.APPLYtr(x,y);
//   cout << y << endl;

//   cout << "- testing non-quadratic lower triangular matrices:" << endl;
//   RawMatrix<double,LowerTriangularArray2D<double> > Lnonq1(3,4,"1 2 3 4 5 6");
//   cout << Lnonq1 << endl;
//   cout << "  (this has " << Lnonq1.size() << " nonzero entries)" << endl;
//   RawMatrix<double,LowerTriangularArray2D<double> > Lnonq2(4,3,"1 2 3 4 5 6 7 8 9");
//   cout << Lnonq2 << endl;
//   cout << "  (this has " << Lnonq2.size() << " nonzero entries)" << endl;

//   cout << "- first matrix, testing APPLY for x=" << x << ", result: ";
//   Lnonq1.APPLY(x,y);
//   cout << y << endl;
//   cout << "- first matrix, testing APPLYtr for x=";
//   x.resize(Lnonq1.rowdim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   Lnonq1.APPLYtr(x,y);
//   cout << y << endl;
//   cout << "- second matrix, testing APPLY for x=";
//   x.resize(Lnonq2.coldim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   Lnonq2.APPLY(x,y);
//   cout << y << endl;
//   cout << "- second matrix, testing APPLYtr for x=";
//   x.resize(Lnonq2.rowdim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   Lnonq2.APPLYtr(x,y);
//   cout << y << endl;

//   cout << "- testing non-quadratic upper triangular matrices:" << endl;
//   RawMatrix<double,UpperTriangularArray2D<double> > Rnonq1(3,4,"1 2 3 4 5 6 7 8 9");
//   cout << Rnonq1 << endl;
//   cout << "  (this has " << Rnonq1.size() << " nonzero entries)" << endl;
//   RawMatrix<double,UpperTriangularArray2D<double> > Rnonq2(4,3,"1 2 3 4 5 6");
//   cout << Rnonq2 << endl;
//   cout << "  (this has " << Rnonq2.size() << " nonzero entries)" << endl;

//   cout << "- first matrix, testing APPLY for x=";
//   x.resize(Rnonq1.coldim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   Rnonq1.APPLY(x,y);
//   cout << y << endl;
//   cout << "- first matrix, testing APPLYtr for x=";
//   x.resize(Rnonq1.rowdim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   Rnonq1.APPLYtr(x,y);
//   cout << y << endl;
//   cout << "- second matrix, testing APPLY for x=";
//   x.resize(Rnonq2.coldim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   Rnonq2.APPLY(x,y);
//   cout << y << endl;
//   cout << "- second matrix, testing APPLYtr for x=";
//   x.resize(Rnonq2.rowdim());
//   for (int i(0); i < x.dim(); i++) x(i)=i+1;
//   cout << x << ", result: ";
//   Rnonq2.APPLYtr(x,y);
//   cout << y << endl;

  return 0;
}

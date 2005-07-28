#include <iostream>
#include "MathTL.h"

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the matrix decomposition routines ..." << endl;

  typedef Matrix<double>::size_type size_type;
  size_type banddim = 8;
  SymmetricMatrix<double> A(banddim);
  LowerTriangularMatrix<double> L;
  for (size_type i(0); i < banddim; i++)
    {
      if (i >= 1) A(i, i-1) = -1;
      A(i, i) = 2;
    }
  cout << "- a (symmetric) band matrix A:" << endl
       << A;

  cout << "- Cholesky decomposition of A ..." << endl;
  bool result = CholeskyDecomposition(A, L);
  if (result)
    {
      cout << "  ... succeeded, L=" << endl << L << endl
	   << "  and L*L^T=" << endl;

      // as long as we don't have matrix multiplication and transposition,
      // we emulate this:
      Matrix<double> M(banddim, banddim);
      for (size_type i(0); i < banddim; i++)
 	for (size_type j(0); j < banddim; j++)
 	  {
 	    for (size_type k(0); k <= std::min(i,j); k++)
	      M(i,j) += L(i,k) * L(j,k);
 	  }
      cout << M;
    }
  else
    cout << "  ... failed!" << endl;
      
  cout << "- perform QU decomposition of A..." << endl;
  QUDecomposition<double> qu(A);
  if (qu.hasFullRank())
    cout << "  * A has full rank!" << endl;
  else
    cout << "  * A does not have full rank!" << endl;
  cout << "  * U:" << endl;
  UpperTriangularMatrix<double> U;
  qu.getU(U);
  cout << U << endl;
  cout << "  * Q:" << endl;
  Matrix<double> Q;
  qu.getQ(Q);
  cout << Q << endl;
  cout << "  * check:" << endl;
  Matrix<double> QU(banddim);
  for (size_type i(0); i < banddim; i++)
    for (size_type j(0); j < banddim; j++)
      {
	for (size_type k(0); k <= j; k++)
	  QU(i,j) += Q(i,k) * U(k,j);
      }
  cout << QU;

  cout << "- perform SVD of A..." << endl;
  SVD<double> svd(A);
  cout << "  * U*S:" << endl;
  Matrix<double> US;
  svd.getUS(US);
  cout << US << endl;
  cout << "  * V:" << endl;
  Matrix<double> V;
  svd.getV(V);
  cout << V << endl;
  cout << "  * check of A-UV*S yields error " << row_sum_norm(Matrix<double>(A)-(US*V)) << endl;
  cout << "  * singular values:" << endl;
  Vector<double> S;
  svd.getS(S);
  cout << S << endl;

  Matrix<double> U2;
  svd.getU(U2);
  for (unsigned int i(0); i < U2.row_dimension(); i++)
    for (unsigned int j(0); j < U2.column_dimension(); j++)
      U2(i, j) *= S[j];
  
  cout << "  * check of A-U*S*V yields error " << row_sum_norm(Matrix<double>(A)-(U2*V)) << endl;

  cout << "- check linear solver using QU decomposition:" << endl
       << "  * a small test matrix B:" << endl;
  Matrix<double> B(3, 2, "1 1 -1 1 0 1");
  cout << B;
  Vector<double> xexact(2);
  xexact[0] = 1; xexact[1] = 2;
  cout << "  * a vector x: " << xexact << endl;
  Vector<double> c(3);
  for (unsigned int i(0); i < 3; i++)
    {
      c[i] = 0;
      for (unsigned int j(0); j < 2; j++)
	c[i] += B(i, j) * xexact[j];
    }
  cout << "  * rhs vector c: " << c << endl;
  Vector<double> x;
  QUDecomposition<double> qu2(B);
  qu2.solve(c, x);
  cout << "  * solve() output x: " << x << endl;

  cout << "- check computation of inverses via QU decomposition:" << endl
       << "  * a small test matrix M:" << endl;
  Matrix<double> M(3, 3, "1 2 3 0 1 2 0 0 1"), MInv;
  cout << M;
  QUDecomposition<double> qu3(M);
  qu3.inverse(MInv);
  cout << "  * inverse() output:" << endl << MInv;

  return 0;
}

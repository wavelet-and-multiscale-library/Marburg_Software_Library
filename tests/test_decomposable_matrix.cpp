#include <iostream>
#include <numerics/decomposable_matrix.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing decomposable matrices ..." << endl;

  DecomposableMatrix<double> M(2);
  cout << "- a default decomposable matrix M:" << endl
       << M;

//   typedef DecomposableMatrix<double>::size_type size_type;
//   size_type banddim = 8;
//   SymmetricMatrix<double> A(banddim);
//   LowerTriangularMatrix<double> L;
//   for (size_type i(0); i < banddim; i++)
//     {
//       if (i >= 1) A(i, i-1) = -1;
//       A(i, i) = 2;
//     }
//   cout << "- a (symmetric) band matrix A:" << endl
//        << A;

//   cout << "- Cholesky decomposition of A ..." << endl;
//   bool result = CholeskyDecomposition(A, L);
//   if (result)
//     {
//       cout << "  ... succeeded, L=" << endl << L << endl
// 	   << "  and L*L^T=" << endl;

//       // as long as we don't have matrix multiplication and transposition,
//       // we emulate this:
//       Matrix<double> M(banddim, banddim);
//       for (size_type i(0); i < banddim; i++)
//  	for (size_type j(0); j < banddim; j++)
//  	  {
//  	    for (size_type k(0); k <= std::min(i,j); k++)
// 	      M(i,j) += L(i,k) * L(j,k);
//  	  }
//       cout << M;
//     }
//   else
//     cout << "  ... failed!" << endl;
      
//   cout << "- perform QR decomposition of A..." << endl;
//   QRDecomposition<double> qr(A);
//   if (qr.hasFullRank())
//     cout << "  * A has full rank!" << endl;
//   else
//     cout << "  * A does not have full rank!" << endl;
//   cout << "  * R:" << endl;
//   UpperTriangularMatrix<double> R;
//   qr.getR(R);
//   cout << R << endl;
//   cout << "  * Q:" << endl;
//   Matrix<double> Q;
//   qr.getQ(Q);
//   cout << Q << endl;
//   cout << "  * check:" << endl;
//   Matrix<double> QR(banddim);
//   for (size_type i(0); i < banddim; i++)
//     for (size_type j(0); j < banddim; j++)
//       {
// 	for (size_type k(0); k <= j; k++)
// 	  QR(i,j) += Q(i,k) * R(k,j);
//       }
//   cout << QR;

//   cout << "- perform SVD of A..." << endl;
//   SVD<double> svd(A);
//   cout << "  * U*S:" << endl;
//   Matrix<double> US;
//   svd.getUS(US);
//   cout << US << endl;
//   cout << "  * V:" << endl;
//   Matrix<double> V;
//   svd.getV(V);
//   cout << V << endl;
//   cout << "  * check of A-UV*S yields error " << row_sum_norm(Matrix<double>(A)-(US*V)) << endl;
//   cout << "  * singular values:" << endl;
//   Vector<double> S;
//   svd.getS(S);
//   cout << S << endl;

//   Matrix<double> U;
//   svd.getU(U);
//   for (unsigned int i(0); i < U.row_dimension(); i++)
//     for (unsigned int j(0); j < U.column_dimension(); j++)
//       U(i, j) *= S[j];
  
//   cout << "  * check of A-U*S*V yields error " << row_sum_norm(Matrix<double>(A)-(U*V)) << endl;

//   cout << "- check linear solver using QR decomposition:" << endl
//        << "  * a small test matrix B:" << endl;
//   Matrix<double> B(3, 2, "1 1 -1 1 0 1");
//   cout << B;
//   Vector<double> xexact(2);
//   xexact[0] = 1; xexact[1] = 2;
//   cout << "  * a vector x: " << xexact << endl;
//   Vector<double> c(3);
//   for (unsigned int i(0); i < 3; i++)
//     {
//       c[i] = 0;
//       for (unsigned int j(0); j < 2; j++)
// 	c[i] += B(i, j) * xexact[j];
//     }
//   cout << "  * rhs vector c: " << c << endl;
//   Vector<double> x;
//   QRDecomposition<double> qr2(B);
//   qr2.solve(c, x);
//   cout << "  * solve() output x: " << x << endl;

//   cout << "- check computation of inverses via QR decomposition:" << endl
//        << "  * a small test matrix M:" << endl;
//   Matrix<double> M(3, 3, "1 2 3 0 1 2 0 0 1"), MInv;
//   cout << M;
//   QRDecomposition<double> qr3(M);
//   qr3.inverse(MInv);
//   cout << "  * inverse() output:" << endl << MInv;

  return 0;
}

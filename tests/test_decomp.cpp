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
      
  cout << "- perform QR decomposition of A..." << endl;
  QRDecomposition<double> qr(A);
  if (qr.hasFullRank())
    cout << "  * A has full rank!" << endl;
  else
    cout << "  * A does not have full rank!" << endl;
  cout << "  * R:" << endl;
  UpperTriangularMatrix<double> R;
  qr.getR(R);
  cout << R << endl;
  cout << "  * Q:" << endl;
  Matrix<double> Q;
  qr.getQ(Q);
  cout << Q << endl;
  cout << "  * check:" << endl;
  Matrix<double> QR(banddim, banddim);
  for (size_type i(0); i < banddim; i++)
    for (size_type j(0); j < banddim; j++)
      {
	for (size_type k(0); k <= j; k++)
	  QR(i,j) += Q(i,k) * R(k,j);
      }
  cout << QR;

//   cout << "- perform SVD of A..." << endl;
//   SVD<double> svd(A);
//   cout << "  * U*S:" << endl;
//   RawMatrix<double> US;
//   svd.getUS(US);
//   cout << US << endl;
//   cout << "  * V:" << endl;
//   RawMatrix<double> V;
//   svd.getV(V);
//   cout << V << endl;
//   cout << "  * check:" << endl
//        << US*transpose(V) << endl;
//   cout << "  * singular values:" << endl;
//   DenseArray1D<double> S;
//   svd.getS(S);
//   cout << S << endl;

  return 0;
}

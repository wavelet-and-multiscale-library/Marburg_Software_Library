#include <iostream>
#include <cmath>
#include "MathTL.h"

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the computation of extremal eigenvalues..." << endl;

//   int banddim = 8;
//   RawMatrix<double,SymmetricArray2D<double> > A(banddim, banddim);
//   for (int i(0); i < banddim; i++)
//     {
//       if (i >= 1) A(i, i-1) = -1;
//       A(i, i) = 2;
//     }
//   cout << "- a (symmetric) band matrix A:" << endl
//        << A << endl;
  
//   RawVector<double> xk, err;
//   cout << "- calculating the maximal eigenvalue of A ..." << endl;
//   int iterations;
//   double lambda = PowerIteration(A, xk, 1e-6, 100, iterations);
//   A.APPLY(xk, err);
//   err -= lambda * xk;
//   cout << "  ... yields (after " << iterations << " iterations) lambda=" << lambda
//        << ", corresponding (1-normalized) eigenvector:" << endl
//        << xk << endl
//        << "  with \\|A*xk-lambdak*xk\\|_\\infty=" << maxnorm(err) << endl;
  
//   cout << "- calculating the minimal eigenvalue of A ..." << endl;
//   lambda = 1./InversePowerIteration(A, xk, 1e-6, 200, iterations);
//   A.APPLY(xk, err);
//   err -= lambda * xk;
//   cout << "  ... yields (after " << iterations << " iterations) lambda=" << lambda
//        << ", corresponding (1-normalized) eigenvector:" << endl
//        << xk << endl
//        << "  with \\|A^{-1}*xk-lambdak*xk\\|_\\infty=" << maxnorm(err) << endl;

//   cout << "- spectral condition number of A (with CondSymm()): "
//        << CondSymm(A) << endl;
//   cout << "- spectral condition number of A (with CondNonSymm()): "
//        << CondNonSymm(A) << endl;

//   cout << "- calculating ALL eigenvalues and eigenvectors of A ..." << endl;
//   RawMatrix<double> evecs;
//   RawVector<double> evals;
//   SymmEigenvalues(A, evals, evecs);
//   cout << "  ... eigenvalues:" << endl << evals << endl
//        << "  ... eigenvectors (row-wise):" << endl << evecs << endl;
  
//   cout << "- check the eigenvectors:" << endl;
//   for (int i(0); i < evals.dim(); i++)
//     {
//       RawVector<double> x(A.coldim());
//       for (int j(0); j < x.dim(); j++)
// 	x(j) = evecs(j,i);
//       A.APPLY(x,err);
//       err -= evals(i)*x;
//       cout << "  vector " << i+1 << " has \\|A*x_i-lambda_i*x_i\\|_\\infty="
// 	   <<  maxnorm(err) << endl;
//     }
  
  return 0;
}

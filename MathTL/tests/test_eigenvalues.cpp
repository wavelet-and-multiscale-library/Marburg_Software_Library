#include <iostream>
#include <cmath>
#include <time.h>
#include "MathTL.h"

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the computation of extremal eigenvalues..." << endl;

  unsigned int banddim = 10;
  SymmetricMatrix<double> A(banddim);
  for (unsigned int i(0); i < banddim; i++)
    {
      if (i >= 1) A(i, i-1) = -1;
      A(i, i) = 2;
    }
  cout << "- a (symmetric) band matrix A:" << endl
       << A;
  
  cout << "- exact eigenvalues of A:" << endl;
  for (unsigned int i(0); i < banddim; i++)
    cout << 2-2*cos((i+1)*M_PI/(banddim+1)) << " "; // siehe Numerik-Uebungen...
  cout << endl;
  
  Vector<double> xk(banddim, false), err(banddim, false);

  srand(time(0)); // !

  for (unsigned int i(0); i < xk.size(); i++) xk[i] = random_double();
//   xk = 1;
  cout << "- calculating the maximal eigenvalue of A ..." << endl;
  unsigned int iterations;
  double lambda = PowerIteration(A, xk, 1e-6, 100, iterations);
  A.apply(xk, err);
  err.add(-lambda, xk);
  cout << "  ... yields (after " << iterations << " iterations) lambda=" << lambda
       << ", corresponding (1-normalized) eigenvector:" << endl
       << xk << endl
       << "  with \\|A*xk-lambdak*xk\\|_\\infty=" << linfty_norm(err) << endl;
  
  cout << "- calculating the minimal eigenvalue of A ..." << endl;
  for (unsigned int i(0); i < xk.size(); i++) xk[i] = random_double();
//   xk = 1;
  lambda = InversePowerIteration(A, xk, 1e-6, 200, iterations);
  A.apply(xk, err);
  err.add(-lambda, xk);
  cout << "  ... yields (after " << iterations << " iterations) lambda=" << lambda
       << ", corresponding (1-normalized) eigenvector:" << endl
       << xk << endl
       << "  with \\|A^{-1}*xk-lambdak*xk\\|_\\infty=" << linfty_norm(err) << endl;

  cout << "- calculating the minimal eigenvalue again, with spectral shift lambda=0.05 ..." << endl;
  for (unsigned int i(0); i < xk.size(); i++) xk[i] = random_double();
//   xk = 1;
  lambda = InversePowerIteration(A, xk, 0.05, 1e-6, 200, iterations);
  A.apply(xk, err);
  err.add(-lambda, xk);
  cout << "  ... yields (after " << iterations << " iterations) lambda=" << lambda
       << ", corresponding (1-normalized) eigenvector:" << endl
       << xk << endl
       << "  with \\|A^{-1}*xk-lambdak*xk\\|_\\infty=" << linfty_norm(err) << endl;  

  cout << "- spectral condition number of A (with CondSymm()): "
       << CondSymm(A) << endl;
  cout << "- spectral condition number of A (with CondNonSymm()): "
       << CondNonSymm(A) << endl;

  cout << "- calculating ALL eigenvalues and eigenvectors of A ..." << endl;
  Matrix<double> evecs;
  Vector<double> evals;
  SymmEigenvalues(A, evals, evecs);
  cout << "  ... eigenvalues:" << endl << evals << endl
       << "  ... eigenvectors (row-wise):" << endl << evecs;
  
  cout << "- check the eigenvectors:" << endl;
  for (unsigned int i(0); i < evals.size(); i++)
    {
      Vector<double> x(A.column_dimension(), false);
      for (unsigned int j(0); j < x.size(); j++)
 	x[j] = evecs(j,i);
      A.apply(x,err);
      err.add(-evals(i), x);
      cout << "  vector " << i+1 << " has \\|A*x_i-lambda_i*x_i\\|_\\infty="
 	   << linfty_norm(err) << endl;
     }

  cout << "- Lanczos iteration for A:" << endl;
  double lambdamin, lambdamax;
  LanczosIteration(A, 1e-5, lambdamin, lambdamax, 100, iterations);
  cout << "  lambdamin=" << lambdamin << ", lambdamax=" << lambdamax
       << ", " << iterations << " iterations needed" << endl;

  return 0;
}

#include <iostream>
#include <cmath>
#include "MathTL.h"

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing the iterative solvers..." << endl;
  
  unsigned int banddim = 10, iterations, maxiter = 2000;
  SymmetricMatrix<double> A(banddim);
  for (unsigned int i(0); i < banddim; i++)
    {
      if (i >= 1) A(i, i-1) = -1;
      A(i, i) = 2;
    }
  cout << "- a (symmetric) band matrix A:" << endl
       << A;
  
  Vector<double> b(banddim, false), xk(banddim, false), err(banddim, false);
  
  for (unsigned int i(0); i < banddim; i++)
    b(i) = i;
  cout << "- the right-hand side b:" << endl
       << b << endl;
  
  for (unsigned int i(0); i <= 3; i++)
    {
      double omega = 0.35 + i*0.05;
      xk = 0;
      cout << "- Richardson iteration with omega=" << omega << " ..." << endl;
      Richardson(A, b, xk, omega, 1e-8, maxiter, iterations);
      A.apply(xk, err);
      err -= b;
      cout << "  ... yields a solution xk=" << endl
	   << xk
	   << "  with \\|A*xk-b\\|_\\infty=" << linfty_norm(err)
	   << " after " << iterations << " iterations." << endl;
    }

//   cout << "- determine optimal relaxation parameter for Richardson iteration:" << endl;
//   double lambdamax = PowerIteration(A, xk, 1e-6, 1000, iterations);
//   double lambdamin = 1./InversePowerIteration(A, xk, 1e-6, 1000, iterations);
//   cout << "  lambdamax=" << lambdamax << ", 2/lambdamax=" << 2./lambdamax
//        << ", lambdamin=" << lambdamin << endl;
//   double omegastar = 2./(lambdamin+lambdamax);
//   xk = 0;
//   cout << "  Richardson iteration with omegastar=" << omegastar << " ..." << endl;
//   Richardson(A, b, xk, omegastar, 1e-8, maxiter, iterations);
//   A.APPLY(xk, err);
//   err -= b;
//   cout << "  ... yields a solution xk=" << endl
//        << xk
//        << "  with \\|A*xk-b\\|_\\infty=" << maxnorm(err)
//        << " after " << iterations << " iterations." << endl;

  xk = 0; xk(0) = 1;
  cout << "- Jacobi iteration ..." << endl;
  Jacobi(A, b, xk, 1e-8, maxiter, iterations);
  A.apply(xk, err);
  err -= b;
  cout << "  ... yields a solution xk=" << endl
       << xk
       << "  with \\|A*xk-b\\|_\\infty=" << linfty_norm(err)
       << " after " << iterations << " iterations." << endl;
  
  xk = 0; xk(0) = 1;
  cout << "- Gauss-Seidel iteration ..." << endl;
  GaussSeidel(A, b, xk, 1e-8, maxiter, iterations);
  A.apply(xk, err);
  err -= b;
  cout << "  ... yields a solution xk=" << endl
       << xk
       << "  with \\|A*xk-b\\|_\\infty=" << linfty_norm(err)
       << " after " << iterations << " iterations." << endl;
  
  xk = 0; xk(0) = 1;
  cout << "- SOR iteration ..." << endl;
  SOR(A, b, xk, 1.2, 1e-8, maxiter, iterations);
  A.apply(xk, err);
  err -= b;
  cout << "  ... yields a solution xk=" << endl
       << xk
       << "  with \\|A*xk-b\\|_\\infty=" << linfty_norm(err)
       << " after " << iterations << " iterations." << endl;
  
  xk = 0; xk(0) = 1;
  cout << "- SSOR iteration ..." << endl;
  SSOR(A, b, xk, 1.2, 1e-8, maxiter, iterations);
  A.apply(xk, err);
  err -= b;
  cout << "  ... yields a solution xk=" << endl
       << xk
       << "  with \\|A*xk-b\\|_\\infty=" << linfty_norm(err)
       << " after " << iterations << " iterations." << endl;

  xk = 0; xk(0) = 1;
  cout << "- CG iteration ..." << endl;
  CG(A, b, xk, 1e-8, maxiter, iterations);
  A.apply(xk, err);
  err -= b;
  cout << "  ... yields a solution xk=" << endl
       << xk << endl
       << "  with \\|A*xk-b\\|_\\infty=" << linfty_norm(err)
       << " after " << iterations << " iterations." << endl;

  return 0;
}

#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>

#include <interval/s_basis.h>
#include <interval/s_support.h>
#include <interval/interval_evaluate.h>
#include <galerkin/gramian.h>
#include <galerkin/galerkin_utils.h>

#include <numerics/eigenvalues.h>
#include <numerics/iteratsolv.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;


/*
  [P] M. Primbs: Stabile biorthogonale Wavelet-Basen auf dem Intervall
      Dissertation, Univ. Duisburg-Essen, 2006
 */


typedef SBasis Basis;


int main()
{
  Basis basis;
  Basis::Index lambda;
  set<Basis::Index> Lambda;

  int jmax = basis.j0()+3;
  for (lambda = basis.first_generator(basis.j0()); lambda <= basis.last_wavelet(jmax); ++lambda)
    Lambda.insert(lambda);
  SparseMatrix<double> G(Lambda.size());

  cout << "Estimating the Riesz constants via Eigenvalues of the Gram Matrix" << endl;
  // cf. [P], Kap. 6
  cout << "- setting up Gramian matrix ..." << endl;
  IntervalGramian<Basis> problem(basis, InfiniteVector<double,Basis::Index>());
  cout << "- solving Eigenvalue problem ..." << endl;
  setup_stiffness_matrix(problem, Lambda, G, true);
  double lambdamin, lambdamax;
  unsigned int iterations;
  LanczosIteration(G, 1e-6, lambdamin, lambdamax, 200, iterations);
  cout << "Needed " << iterations << " iterations for the Eigenvalue problem" << endl;
  cout << "C_1 = " << sqrt(lambdamin) << endl;
  cout << "C_2 = " << sqrt(lambdamax) << endl;
  cout << "kappa = C_2^2/C_1^2 = " << lambdamax/lambdamin << endl;

  return 0;
}

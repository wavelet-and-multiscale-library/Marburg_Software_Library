// choose which basis to use
#define BASIS_S 0
#define BASIS_P 1
#define BASIS_DS 2
#define BASIS_JL 3

#define BASIS BASIS_S


#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>

#if (BASIS == BASIS_S)
 #include <interval/s_basis.h>
 #include <interval/s_support.h>
 #include <interval/interval_evaluate.h>
 typedef WaveletTL::SBasis Basis;
#elif (BASIS == BASIS_P)
 #include <interval/p_basis.h>
 #include <interval/p_support.h>
 #include <interval/p_evaluate.h>
 typedef WaveletTL::PBasis<4, 4> Basis;
#elif (BASIS == BASIS_DS)
 #include <interval/ds_basis.h>
 #include <interval/ds_support.h>
 #include <interval/ds_evaluate.h>
 typedef WaveletTL::DSBasis<4, 4> Basis;
#elif (BASIS == BASIS_JL)
 #include <interval/jl_basis.h>
 #include <interval/jl_support.h>
 #include <interval/jl_evaluate.h>
 typedef WaveletTL::JLBasis Basis;
#endif

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


int main()
{
  #if ( (BASIS == BASIS_P) || (BASIS == BASIS_DS) )
  Basis basis(2, 2);
  #else
  Basis basis;
  #endif
  Basis::Index lambda;
  set<Basis::Index> Lambda;

  int jmax = 12; //basis.j0()+3;
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

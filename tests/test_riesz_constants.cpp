// choose which basis to use
#define BASIS_S 0
#define BASIS_P 1
#define BASIS_DS 2
#define BASIS_JL 3

#define BASIS BASIS_S

#define JMAX_START 18 // basis.j0()
#define JMAX_END 25 // 16
#define JMAX_STEP 1

#define SAVE_RESULTS


#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>

#if (BASIS == BASIS_S)
 #include <interval/s_basis.h>
 #include <interval/s_support.h>
 #include <interval/interval_evaluate.h>
 typedef WaveletTL::SBasis Basis;
 #define BASIS_NAME "s"
#elif (BASIS == BASIS_P)
 #include <interval/p_basis.h>
 #include <interval/p_support.h>
 #include <interval/p_evaluate.h>
 typedef WaveletTL::PBasis<4, 4> Basis;
 #define BASIS_NAME "p"
#elif (BASIS == BASIS_DS)
 #include <interval/ds_basis.h>
 #include <interval/ds_support.h>
 #include <interval/ds_evaluate.h>
 typedef WaveletTL::DSBasis<4, 4> Basis;
 #define BASIS_NAME "ds"
#elif (BASIS == BASIS_JL)
 #include <interval/jl_basis.h>
 #include <interval/jl_support.h>
 #include <interval/jl_evaluate.h>
 typedef WaveletTL::JLBasis Basis;
 #define BASIS_NAME "jl"
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

  cout << "Estimating the Riesz constants via Eigenvalues of the Gramian matrix" << endl;
  // cf. [P], Kap. 6
  #ifdef SAVE_RESULTS
  std::ofstream fs;
  ostringstream filename;
  filename << "riesz_constants_" << BASIS_NAME << ".dat";
  fs.open(filename.str().c_str());
  fs << "BASIS = " << BASIS_NAME << endl;
  fs << "jmax\tC_1\tC_2\tkappa" << endl;
  #endif // SAVE_RESULTS
  for (int jmax = JMAX_START; jmax <= JMAX_END; jmax += JMAX_STEP) {
    cout << "jmax = " << jmax << endl;
    Lambda.clear();
    for (lambda = basis.first_generator(basis.j0()); lambda <= basis.last_wavelet(jmax); ++lambda)
      Lambda.insert(lambda);
    SparseMatrix<double> G(Lambda.size());

    cout << "- setting up Gramian matrix ..." << endl;
    IntervalGramian<Basis> problem(basis, InfiniteVector<double,Basis::Index>());
    setup_stiffness_matrix(problem, Lambda, G, true);
    cout << "- preconditioning Gramian matrix ..." << endl;
    Vector<double> diag(G.row_dimension());
    for (unsigned int i = 0; i < G.row_dimension(); i++)
      diag[i] = G.get_entry(i, i);
    SparseMatrix<double> H(G.row_dimension(), G.column_dimension());
    for (unsigned int row = 0; row < G.row_dimension(); row++) {
      for (unsigned int col = 0; col < G.column_dimension(); col++) {
        double entry = G.get_entry(row, col);
        if (entry != 0) {
          H.set_entry(row, col, entry/(sqrt(diag[row])*sqrt(diag[col])));
        }
      }
    }
    cout << "- solving Eigenvalue problem ..." << endl;
    double lambdamin, lambdamax;
    unsigned int iterations;
    LanczosIteration(H, 1e-10, lambdamin, lambdamax, 200, iterations);
    cout << "Needed " << iterations << " iterations for the Eigenvalue problem" << endl;
    cout << "C_1 = " << sqrt(lambdamin) << endl;
    cout << "C_2 = " << sqrt(lambdamax) << endl;
    cout << "kappa = C_2^2/C_1^2 = " << lambdamax/lambdamin << endl;
    #ifdef SAVE_RESULTS
    fs << jmax << "\t" << sqrt(lambdamin) << "\t" << sqrt(lambdamax) << "\t" << lambdamax/lambdamin << endl;
    #endif // SAVE_RESULTS
  }
  #ifdef SAVE_RESULTS
  fs.close();
  #endif // SAVE_RESULTS

  return 0;
}

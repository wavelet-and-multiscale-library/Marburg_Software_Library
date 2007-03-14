#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/vector_norms.h>
#include <algebra/matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix_norms.h>
#include <utils/random.h>
#include <ctime>

#include <interval/s_basis.h>
#include <interval/s_support.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;


int main()
{
  cout << "Testing the S bases ..." << endl << endl;

  SBasis basis;

  cout << "first generator of level j_0: " << basis.first_generator(basis.j0()) << endl;
  cout << "last generator of level j_0: " << basis.last_generator(basis.j0()) << endl;
  cout << "first wavelet of level j_0: " << basis.first_wavelet(basis.j0()) << endl;
  cout << "last wavelet of level j_0: " << basis.last_wavelet(basis.j0()) << endl;
  cout << basis.DeltaRmax(basis.j0()+1) << endl;

  Matrix<double> Mj0 = Matrix<double>(30,14);
  Matrix<double> MTj0 = Matrix<double>(30,14);
  Matrix<double> Mj1 = Matrix<double>(30,16);
  Matrix<double> MTj1 = Matrix<double>(30,16);
  SparseMatrix<double> Gj = SparseMatrix<double>(30,30);
  Matrix<double> GTj = Matrix<double>(30,30);
  InfiniteVector<double, SBasis::Index> coeff;
  SBasis::Index lambda;
  unsigned int row, col;

  cout << endl << "Constructing M_{j,0} via refine_1" << endl;
  for (lambda = basis.first_generator(basis.j0()); lambda <= basis.last_generator(basis.j0()); ++lambda) {
    basis.refine_1(lambda, coeff);
    col = (lambda.k()-basis.DeltaLmin())*2+lambda.c();
    for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it) {
      Mj0.set_entry((it.index().k()-basis.DeltaLmin())*2+it.index().c(), col, *it);
    }
  }
  cout << Mj0 << endl;

  cout << "Constructing \\tilde M_{j,0} via refine_t_1" << endl;
  for (lambda = basis.first_generator(basis.j0()); lambda <= basis.last_generator(basis.j0()); ++lambda) {
    basis.refine_t_1(lambda, coeff);
    col = (lambda.k()-basis.DeltaLTmin())*2+lambda.c();
    for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it) {
      MTj0.set_entry((it.index().k()-basis.DeltaLTmin())*2+it.index().c(), col, *it);
    }
  }
  cout << MTj0 << endl;

  cout << "Constructing M_{j,1} via reconstruct_1" << endl;
  for (lambda = basis.first_wavelet(basis.j0()); lambda <= basis.last_wavelet(basis.j0()); ++lambda) {
    basis.reconstruct_1(lambda, coeff);
    col = (lambda.k()-basis.Nablamin())*2+lambda.c();
    for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it) {
//      cout << it.index() << ", " << (it.index().k()-basis.DeltaLmin())*2+it.index().c() << " ";
      Mj1.set_entry((it.index().k()-basis.DeltaLmin())*2+it.index().c(), col, *it);
    }
  }
  cout << Mj1 << endl;

  cout << "Constructing \\tilde M_{j,1} via reconstruct_t_1" << endl;
  for (lambda = basis.first_wavelet(basis.j0()); lambda <= basis.last_wavelet(basis.j0()); ++lambda) {
    basis.reconstruct_t_1(lambda, coeff);
    col = (lambda.k()-basis.Nablamin())*2+lambda.c();
    for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it) {
      MTj1.set_entry((it.index().k()-basis.DeltaLTmin())*2+it.index().c(), col, *it);
    }
  }
  cout << MTj1 << endl;

  cout << "Constructing G_j ( = (\\tilde M_{j,0}, \\tilde M_{j,1})^T ) via decompose_1" << endl;
  for (lambda = basis.first_generator(basis.j0()+1); lambda <= basis.last_generator(basis.j0()+1); ++lambda) {
    basis.decompose_1(lambda, coeff);
    col = (lambda.k()-basis.DeltaLmin())*2+lambda.c();
    for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it) {
      if (it.index().e() == E_GENERATOR) // generator part
        row = (it.index().k() - basis.DeltaLmin())*2 + it.index().c();
      else // wavelet part
        row = 2*basis.DeltaRmax(basis.j0()) + (it.index().k() - basis.Nablamin())*2 + it.index().c();
      Gj.set_entry(row, col, *it);
    }
  }
  cout << Gj << endl;

  cout << "Constructing \\tilde G_j ( = (M_{j,0}, M_{j,1})^T ) via decompose_t_1" << endl;
  for (lambda = basis.first_generator(basis.j0()+1); lambda <= basis.last_generator(basis.j0()+1); ++lambda) {
//    cout << "lambda = " << lambda << endl;
    basis.decompose_t_1(lambda, coeff);
    col = (lambda.k()-basis.DeltaLTmin())*2+lambda.c();
    for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it) {
      if (it.index().e() == E_GENERATOR) // generator part
        row = (it.index().k() - basis.DeltaLTmin())*2 + it.index().c();
      else // wavelet part
        row = 2*basis.DeltaRTmax(basis.j0()) + (it.index().k() - basis.Nablamin())*2 + it.index().c();
//      cout << it.index() << " -> (" << row << ", " << col << ")" << endl;
      GTj.set_entry(row, col, *it);
    }
  }
  cout << GTj << endl;


  int j;
  SparseMatrix<double> Mj;
  cout << "Checking M_j . G_j = I" << endl;
  for (j = basis.j0(); j <= basis.j0()+3; j++) {
    Mj = SparseMatrix<double>(basis.Deltasize(j+1)*basis.number_of_components);
    Gj = SparseMatrix<double>(basis.Deltasize(j+1)*basis.number_of_components);
    // setup M_j via refine_1, reconstruct_1
    for (lambda = basis.first_generator(j); lambda <= basis.last_generator(j); ++lambda) {
      basis.refine_1(lambda, coeff);
      col = (lambda.k()-basis.DeltaLmin())*2+lambda.c();
      for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it)
        Mj.set_entry((it.index().k()-basis.DeltaLmin())*2+it.index().c(), col, *it);
    }
    for (lambda = basis.first_wavelet(j); lambda <= basis.last_wavelet(j); ++lambda) {
      basis.reconstruct_1(lambda, coeff);
      col = basis.Deltasize(j)*basis.number_of_components + (lambda.k()-basis.Nablamin())*2+lambda.c();
      for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it)
        Mj.set_entry((it.index().k()-basis.DeltaLmin())*2+it.index().c(), col, *it);
    }
    // setup G_j = (\tilde M_{j,0}, \tilde M_{j,1})^T via decompose_1
    for (lambda = basis.first_generator(j+1); lambda <= basis.last_generator(j+1); ++lambda) {
      basis.decompose_1(lambda, coeff);
      col = (lambda.k()-basis.DeltaLTmin())*2+lambda.c();
      for (InfiniteVector<double, SBasis::Index>::const_iterator it(coeff.begin()), itend(coeff.end()); it != itend; ++it) {
        if (it.index().e() == E_GENERATOR) // generator part
          row = (it.index().k() - basis.DeltaLTmin())*2 + it.index().c();
        else // wavelet part
          row = basis.DeltaRmax(j)*basis.number_of_components + (it.index().k() - basis.Nablamin())*2 + it.index().c();
        Gj.set_entry(row, col, *it);
      }
    }
    // calculate M_j . G_j - I
    SparseMatrix<double> A(Mj*Gj);
    for (row = 0; row < A.row_dimension(); row++) // subtract identity matrix
      A.set_entry(row,row, A.get_entry(row,row) - 1.0);
    // compute norm of difference
    cout << "j = " << j << ": \\|M_j G_j - I\\| = " << maximum_norm(A) << endl;
//    cout << Mj << endl << Gj << endl;
  }
  

  InfiniteVector<double, SBasis::Index> coeff_orig, coeff_decompose, coeff_reconstruct, coeff_diff;
  cout << "decomposing a random coefficient ..." << endl;
  // initialize random generator
  #ifdef HAVE_CLOCK_GETTIME
  struct timespec d;
  clock_gettime(CLOCK_REALTIME, &d);
  srand(d.tv_nsec);
  #else // !HAVE_CLOCK_GETTIME
  time_t rawtime;
  tm* gt;
  rawtime = time(NULL);
  gt = gmtime(&rawtime);
  srand(gt->tm_sec);
  #endif // !HAVE_CLOCK_GETTIME
  // decompose and reconstruct a random coefficient on several levels
  for (j = basis.j0() + 1; j < basis.j0() + 4; j++) {
    lambda = SBasis::Index(j, E_GENERATOR, random_integer(basis.DeltaLmin(),basis.DeltaRmax(j)), random_integer(0, basis.number_of_components-1), &basis);
    cout << "lambda = " << lambda << endl;
//    cout << "after decomposition to j_0:" << endl;
    basis.decompose_1(lambda, basis.j0(), coeff_decompose);
//    cout << coeff_decompose << endl;
    cout << "and after reconstruction:" << endl;
    basis.reconstruct(coeff_decompose, j, coeff_reconstruct);
    coeff_reconstruct.compress();
    cout << coeff_reconstruct << endl;
  }

  cout << "decomposing a random coefficient set ..." << endl;
  // setup random coefficient set on level j
  for (j = 4; j <= 7; j++) {
    cout << "random coefficient set on level j = " << j << ":" << endl;
    coeff_orig.clear();
    for (lambda = basis.first_generator(j); lambda <= basis.last_generator(j); ++lambda)
      coeff_orig.set_coefficient(lambda, random_double());
//  cout << coeff_orig << endl;
    // decompose
//    cout << "after decomposition to j_0:" << endl;
    basis.decompose(coeff_orig, basis.j0(), coeff_decompose);
//    cout << coeff_decompose << endl;
    // reconstruct
//    cout << "and after reconstruction:" << endl;
    basis.reconstruct(coeff_decompose, j, coeff_reconstruct);
//    cout << coeff_reconstruct << endl;
    cout << "Difference of the reconstructed coefficient set to the original one: ";
    coeff_diff = coeff_orig - coeff_reconstruct;
    cout << linfty_norm(coeff_diff) << endl;
  }

  cout << endl << "Testing the (class-member) support routines ..." << endl;
  int k1, k2;
  j = basis.j0()+1;
  lambda = SBasis::Index(j, E_GENERATOR, random_integer(basis.DeltaLmin(),basis.DeltaRmax(j)), random_integer(0, basis.number_of_components-1), &basis); // a random generator index at level j
  basis.primal_support(lambda, k1, k2);
  cout << "support of \\Phi_" << lambda << ": 2^{-" << lambda.j() << "} [" << k1 << ", " << k2 << "]" << endl;
  basis.dual_support(lambda, k1, k2);
  cout << "support of \\tilde\\Phi_" << lambda << ": 2^{-" << lambda.j() << "} [" << k1 << ", " << k2 << "]" << endl;
  lambda = SBasis::Index(j, E_WAVELET, random_integer(basis.Nablamin(),basis.Nablamax(j)), random_integer(0, basis.number_of_components-1), &basis); // a random wavelet index at level j
  basis.primal_support(lambda, k1, k2);
  cout << "support of \\Psi_" << lambda << ": 2^{-" << lambda.j()+lambda.e() << "} [" << k1 << ", " << k2 << "]" << endl;
  basis.dual_support(lambda, k1, k2);
  cout << "support of \\tilde\\Psi_" << lambda << ": 2^{-" << lambda.j()+lambda.e() << "} [" << k1 << ", " << k2 << "]" << endl;

  cout << "support of primal boundary wavelets: ";
  j = basis.j0();
  lambda = SBasis::Index(j, E_WAVELET, basis.Nablamin(), 0, &basis);
  basis.primal_support(lambda, k1, k2);
  cout << lambda << ": 2^{-" << lambda.j()+lambda.e() << "} [" << k1 << ", " << k2 << "]; ";
  lambda = SBasis::Index(j, E_WAVELET, basis.Nablamax(j), 0, &basis);
  basis.primal_support(lambda, k1, k2);
  cout << lambda << ": 2^{-" << lambda.j()+lambda.e() << "} [" << k1 << ", " << k2 << "]" << endl;
  cout << "support of dual boundary generators: ";
  lambda = SBasis::Index(j, E_GENERATOR, basis.DeltaLTmin(), 0, &basis);
  basis.dual_support(lambda, k1, k2);
  cout << lambda << ": 2^{-" << lambda.j() << "} [" << k1 << ", " << k2 << "]; ";
  lambda = SBasis::Index(j, E_GENERATOR, basis.DeltaRTmax(j), 0, &basis);
  basis.dual_support(lambda, k1, k2);
  cout << lambda << ": 2^{-" << lambda.j() << "} [" << k1 << ", " << k2 << "]" << endl;

  return 0;
}

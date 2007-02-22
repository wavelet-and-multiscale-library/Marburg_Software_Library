#include <iostream>
#include <fstream>

#include <algebra/infinite_vector.h>
#include <algebra/matrix.h>

#include <interval/s_basis.h>

using namespace std;
using namespace WaveletTL;

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
  Matrix<double> Gj = Matrix<double>(30,30);
  Matrix<double> GTj = Matrix<double>(30,30);
  InfiniteVector<double, SBasis::Index> coeff;
  SBasis::Index lambda;
  int row, col;

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
    cout << "lambda = " << lambda << endl;
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

  return 0;
}

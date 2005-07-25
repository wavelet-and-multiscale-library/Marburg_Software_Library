#include <iostream>

#include <algebra/infinite_vector.h>
#include <utils/array1d.h>

#include <interval/dku_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DKU bases..." << endl;

  const int d = 2;
  const int dT = 4;
  typedef DKUBasis<d, dT> Basis;
  Basis basis; // Bernstein SVD
//   Basis basis(none);
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;

  cout << "- the default wavelet index: " << Basis::Index(&basis) << endl;

  cout << "- leftmost generator on the coarsest level: " << basis.firstGenerator(basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << basis.lastGenerator(basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << basis.firstWavelet(basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << basis.lastWavelet(basis.j0()) << endl;

#if 0
  cout << "- iterating from first generator on coarsest level to last wavelet on next level:" << endl;
  Basis::Index index(basis.firstGenerator(basis.j0()));
  for (;; ++index) {
    cout << index << endl;
    if (index == basis.lastWavelet(basis.j0()+1)) break;
  }
#endif

  cout << "- index of leftmost right boundary generator on scale j0: "
       << basis.DeltaRmin(basis.j0()) << endl;

#if 0
  cout << "- evaluating some primal generators:" << endl;
  Basis::Index lambda(basis.firstGenerator(basis.j0()));
  for (;; ++lambda) {
    cout << lambda << endl;
    basis.evaluate(lambda, true, 5).matlab_output(cout);
    if (lambda == basis.lastGenerator(basis.j0())) break;
  }
#endif

#if 0
  cout << "- evaluating some dual generators:" << endl;
  lambda = basis.firstGenerator(basis.j0());
  for (;; ++lambda) {
    basis.evaluate(lambda, false, 6).matlab_output(cout);
    if (lambda == basis.lastGenerator(basis.j0())) break;
  }
#endif

  cout << "* another basis:" << endl;

  const int d2 = 3;
  const int dT2 = 5;
  DKUBasis<d2, dT2> basis2;

  cout << "- the (" << d2 << "," << dT2 << ") basis has j0=" << basis2.j0() << endl;

  cout << "- the default wavelet index: " << DKUBasis<d2, dT2>::Index(&basis2) << endl;

  cout << "- leftmost generator on the coarsest level: " << basis2.firstGenerator(basis2.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << basis2.lastGenerator(basis2.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << basis2.firstWavelet(basis2.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << basis2.lastWavelet(basis2.j0()) << endl;

  return 0;
}

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
  const int dT = 2;
  typedef DKUBasis<d, dT> Basis;
  Basis basis(none);
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;

  cout << "- the default wavelet index: " << Basis::Index(&basis) << endl;

  cout << "- leftmost generator on the coarsest level: " << basis.firstGenerator(basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << basis.lastGenerator(basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << basis.firstWavelet(basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << basis.lastWavelet(basis.j0()) << endl;

  cout << "- iterating from first generator on coarsest level to last wavelet on next level:" << endl;
  Basis::Index index(basis.firstGenerator(basis.j0()));
  for (;; ++index) {
    cout << index << endl;
    if (index == basis.lastWavelet(basis.j0()+1)) break;
  }

  return 0;
}

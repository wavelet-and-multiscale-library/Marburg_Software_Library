#include <iostream>
#include <interval/dku.h>
#include <algebra/infinite_vector.h>
#include <utils/array1d.h>

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

//   InfiniteVector<double, Basis::Index> coeff;
//   coeff[RIndex(2,0,0)] = 1.0;
//   coeff[RIndex(2,0,1)] = -3.14;
//   coeff[RIndex(2,0,3)] = 2.0;
//   cout << "  * a small but nontrivial coefficient set:" << endl << coeff;

  return 0;
}

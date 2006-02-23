#include <iostream>
#include <interval/p_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases from [P] ..." << endl;

  const int d = 2;
  const int dT = 4;

  typedef PBasis<d,dT> Basis;
  typedef Basis::Index Index;

  Basis basis(0, 0, 0, 0); // no b.c.'s

  cout << "- d=" << d << ", dT=" << dT << endl;
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;

  return 0;
}

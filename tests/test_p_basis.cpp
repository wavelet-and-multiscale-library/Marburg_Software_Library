#include <iostream>
#include <interval/p_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases from [P] ..." << endl;

  const int d = 2;
  const int dT = 2;

  typedef PBasis<d,dT> Basis;
//   typedef Basis::Index Index;

  Basis basis(0, 0, 0, 0); // no b.c.'s

  return 0;
}

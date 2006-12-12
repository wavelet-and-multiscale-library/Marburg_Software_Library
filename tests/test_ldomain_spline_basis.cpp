#include <iostream>

#include <Ldomain/ldomain_spline_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the L-shaped domain..." << endl;

  const int d  = 2;
  const int dT = 2;

  typedef LDomainSplineBasis<d,dT> Basis;
  Basis basis;

  return 0;
}

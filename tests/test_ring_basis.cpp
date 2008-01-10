#include <iostream>
#include <interval/periodic.h>
#include <ring/ring_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;

  typedef RingBasis<d,dt> Basis;
  Basis basis;

//   typedef Basis::Index Index;

  cout << "- j0=" << basis.j0() << endl;

  return 0;
}

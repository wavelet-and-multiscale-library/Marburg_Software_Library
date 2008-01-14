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

  typedef RingBasis<d,dt,1,0> Basis;
  Basis basis;

  typedef Basis::Index Index;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- first generator on the coarsest level: " << basis.first_generator(basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << basis.last_generator(basis.j0()) << endl;
  cout << "- first wavelet on the coarsest level: " << basis.first_wavelet(basis.j0()) << endl;
  cout << "- last wavelet on the coarsest level: " << basis.last_wavelet(basis.j0()) << endl;

  return 0;
}

#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>

#include <interval/ds_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DS bases..." << endl;

  const int d = 2;
  const int dT = 2;

  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;

//   Basis basis;
  Basis basis(One);
  
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator<Basis>(basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator<Basis>(basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet<Basis>(basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet<Basis>(basis.j0()) << endl;
}

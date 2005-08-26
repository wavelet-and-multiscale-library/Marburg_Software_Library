#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>

#include <Rd/cdf_utils.h>
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

//   Basis basis; // Z={0,1}
  Basis basis(1, 0, 0, 1); // Z={0}
//   Basis basis(0, 1, 1, 0); // Z={1}
//   Basis basis(0, 0, 1, 1); // Z={}
//   Basis basis(0, 0, 0, 0); // should work, DKU basis without b.c.'s at all
  
  cout << "- d=" << d << ", dT=" << dT << endl;
  cout << "- ell1=" << ell1<d>() << ", ell2=" << ell2<d>()
       << ", ell1T=" << ell1T<d,dT>() << ", ell2T=" << ell2T<d,dT>() << endl;
  cout << "- ellT_l=" << basis.ellT_l() << ", ellT_r=" << basis.ellT_r()
       << ", ell_l=" << basis.ell_l() << ", ell_r=" << basis.ell_r() << endl;

  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;
}

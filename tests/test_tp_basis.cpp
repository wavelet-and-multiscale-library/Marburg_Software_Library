#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>

#include <interval/ds_basis.h>
#include <generic/tp_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing tensor product wavelet bases..." << endl;

  typedef DSBasis<2,2> Basis1;
  typedef Basis1::Index Index1;

  typedef DSBasis<2,2> Basis2;
  typedef Basis2::Index Index2;

  typedef TensorProductBasis<Basis1,Basis2> Basis;
  typedef Basis::Index Index;
  Basis basis;
  
  cout << "* a tensor product of 2 DSBasis<2,2> bases:" << endl;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- first generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
  cout << "- first wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
  cout << "- last wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;


}

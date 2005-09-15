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
  cout << "Testing tensor product wavelet bases bases..." << endl;

  typedef DSBasis<2,2> Basis1;
  typedef Basis1::Index Index1;

  typedef DSBasis<2,2> Basis2;
  typedef Basis2::Index Index2;

  TPBasis<Basis1,Basis2> tbasis;
}

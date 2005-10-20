#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>
#include <utils/fixed_array1d.h>

#include <interval/ds_basis.h>
#include <cube/cube_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the cube..." << endl;

  typedef DSBasis<2,2> Basis1D;
  typedef CubeBasis<Basis1D> Basis;
  typedef Basis::Index Index;
  

  Basis basis;
  CubeIndex<Basis1D,2> index(&basis);
  for (int i = 0; i < 255; i++) 
    {
      cout << index << endl;
      ++index;
    }
}

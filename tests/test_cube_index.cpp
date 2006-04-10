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

#if 0
  typedef CubeBasis<Basis1D,1> Basis;
  
  //primal boundary conditions
  FixedArray1D<int,2*1> bound_1_1D;
  bound_1_1D[0] = 1;
  bound_1_1D[1] = 1;

  //dual boundary conditions
  FixedArray1D<int,2*1> bound_1T_1D;
  bound_1T_1D[0] = 0;
  bound_1T_1D[1] = 0;

  Basis basis(bound_1_1D, bound_1T_1D);

  CubeIndex<Basis1D,1,Basis> index(&basis);
  for (int i = 0; i < 64; i++) 
    {
      cout << index << endl;
      ++index;
    }
#endif

  typedef CubeBasis<Basis1D,2> Basis2;
  Basis2 basis2;
  CubeIndex<Basis1D,2> index2(&basis2);
  for (int i = 0; i < 255; i++) {
      cout << index2 << endl;
      ++index2;
    }
  
//   for (Index lambda(first_generator<Basis1D,2,Basis>(&basis, basis.j0()));; ++lambda) {
//     cout << lambda << endl;
//     if (lambda == last_wavelet<Basis1D,2,Basis>(&basis, basis.j0()+1)) break;
//   }

}

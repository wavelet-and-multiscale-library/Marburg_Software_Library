#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>
#include <utils/fixed_array1d.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/jl_basis.h>
#include <interval/jl_support.h>
#include <interval/jl_evaluate.h>

#include <interval/i_index.h>
#include <cube/cube_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the index class of wavelet bases on the cube..." << endl;

//   typedef DSBasis<3,3> Basis1D;
  typedef PBasis<3,3> Basis1D;
//   typedef JLBasis Basis1D; // does not work so far!

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

  //Basis basis(bound_1_1D, bound_1T_1D);
  Basis basis(bound_1_1D);
  //Basis basis;
 
  CubeIndex<Basis1D,1,Basis> index(&basis);
  cout << "* the first 64 indices:" << endl;
  for (int i = 0; i < 64; i++) 
    {
      cout << index << endl;
      ++index;
    }

  typedef Basis1D::Index Index1D;
  
  Basis1D basis1D(bound_1_1D[0], bound_1_1D[1]);

  unsigned int i = 0;
  for (Index1D lambda(first_generator<Basis1D>(&basis1D, basis1D.j0()));; ++lambda) {
    cout << "-----------------------------------------" << endl;
    lambda.set_number();
    cout << lambda << " number  = " << lambda.number() << endl;
    cout << " index  = " << IntervalIndex<Basis1D>(i, &basis1D) << endl;
    if (! (Index1D(i, &basis1D) == lambda))
      abort();
    i++;
    if (lambda == last_wavelet<Basis1D>(&basis1D, basis1D.j0()+2)) break;
  }


  for (int j = basis1D.j0(); j < basis1D.j0()+5; j++ ) {
    if (j == basis1D.j0()) {
      cout << "number of first generator on coarsest level = " << first_generator_num<Basis1D>(&basis1D) << endl;
      cout << "number of last generator on coarsest level = " << last_generator_num<Basis1D>(&basis1D) << endl;
    }
    
    cout << "number of first wavelet on level " << j << " = " << first_wavelet_num<Basis1D>(&basis1D, j) << endl;
    cout << "number of last wavelet on level " << j << " = " << last_wavelet_num<Basis1D>(&basis1D, j) << endl;   
    
  }

//   unsigned int i = 0;
//   for (CubeIndex<Basis1D,1> lambda(first_generator<Basis1D,1,Basis>(&basis, basis.j0()));; ++lambda) {
//     cout << "-----------------------------------------" << endl;
//     cout << lambda << " number  = " << lambda.number() << endl;
//     cout << " index  = " << CubeIndex<Basis1D,1>(i, &basis) << endl;
//      if (! (CubeIndex<Basis1D,1>(i, &basis) == lambda))
//       abort();
//     i++;
//     if (lambda == last_wavelet<Basis1D,1,Basis>(&basis, basis.j0()+2)) break;
//   }

#else

  //primal boundary conditions
  FixedArray1D<int,2*2> bound_1_2D;
  bound_1_2D[0] = 1;
  bound_1_2D[1] = 0;
  bound_1_2D[2] = 0;
  bound_1_2D[3] = 1;
  cout << "testing inhomogeneous boundary conditions bc = " << bound_1_2D << endl;

  //dual boundary conditions
  FixedArray1D<int,2*2> bound_1T_2D;
  bound_1T_2D[0] = 0;
  bound_1T_2D[1] = 0;
  bound_1T_2D[2] = 0;
  bound_1T_2D[3] = 0;

  typedef CubeBasis<Basis1D,2> Basis2;
  //Basis2 basis2(bound_1_2D, bound_1T_2D);
  Basis2 basis2(bound_1_2D);
  //Basis2 basis2;

//   CubeIndex<Basis1D,2> index2(&basis2);
//   for (int i = 0; i < 255; i++) {
//       cout << index2 << endl;
//       ++index2;
//   }

  int i = 0;
  for (CubeIndex<Basis1D,2> lambda(first_generator<Basis1D,2,Basis2>(&basis2, basis2.j0()));; ++lambda) {
    lambda.set_number();
    cout << " index  = " << CubeIndex<Basis1D,2>(i, &basis2);
    cout << " lambda = " << lambda << " number  = " << lambda.number();
    if ((CubeIndex<Basis1D,2>(i, &basis2) != lambda) || (lambda.number()!=i)) cout << " !!! Error !!! ";
    cout << endl;
    if (! (CubeIndex<Basis1D,2>(i, &basis2) == lambda))
      abort();
    i++;
    if (lambda == last_wavelet<Basis1D,2,Basis2>(&basis2, basis2.j0()+1)) break;
    if (i>200) break; 
  }
  
  for (int j = basis2.j0(); j < basis2.j0()+5; j++ ) {
    if (j == basis2.j0()) {
      cout << "number of first generator on coarsest level = " << first_generator_num<Basis1D,2,Basis2>(&basis2) << endl;
      cout << "number of last generator on coarsest level = " << last_generator_num<Basis1D,2,Basis2>(&basis2) << endl;
    }
    
    cout << "number of first wavelet on level " << j << " = " << first_wavelet_num<Basis1D,2,Basis2>(&basis2, j) << endl;
    cout << "number of last wavelet on level " << j << " = " << last_wavelet_num<Basis1D,2,Basis2>(&basis2, j) << endl;   
    
  }

#endif
}

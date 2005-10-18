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
//   Basis basis;
  FixedArray1D<int,4> s, sT;
  s[0] = 1;
  s[1] = 1;
  s[2] = 1;
  s[3] = 1;
  sT[0] = 0;
  sT[1] = 0;
  sT[2] = 0;
  sT[3] = 0;
  Basis basis(s, sT);

  cout << "* a 2D cube basis of DSBasis<2,2> bases:" << endl;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- the default wavelet index w.r.t. the cube basis: " << Index(&basis) << endl;
  cout << "- first generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
  cout << "- first wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
  cout << "- last wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;

#if 0
  cout << "- testing iterator functionality:" << endl;
  for (Index index(first_generator(&basis, basis.j0()));; ++index) {
    cout << index << endl;
    if (index == last_wavelet(&basis, basis.j0()+1)) break;
  }
#endif

#if 1
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(first_generator(&basis, level));
      for (;; ++index)
	{
 	  InfiniteVector<double, Index> origcoeff;
 	  origcoeff[index] = 1.0;
	  
//  	  cout << "* original coeffs:" << endl << origcoeff;

 	  InfiniteVector<double, Index> wcoeff;
 	  basis.decompose(origcoeff, basis.j0(), wcoeff);
	  
//  	  cout << "* after decompose():" << endl << wcoeff;

 	  InfiniteVector<double, Index> transformcoeff;
 	  basis.reconstruct(wcoeff, level, transformcoeff);
	  
// 	  cout << "* after reconstruct():" << endl << transformcoeff;

 	  cout << "* generator: " << index
 	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == last_generator(&basis, level)) break;
	}
    }
#endif
  
#if 0
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose_t() and reconstruct_t() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(first_generator(&basis, level));
      for (;; ++index)
	{
	  InfiniteVector<double, Index> origcoeff;
	  origcoeff[index] = 1.0;
	  
	  InfiniteVector<double, Index> wcoeff;
	  basis.decompose_t(origcoeff, basis.j0(), wcoeff);
	  
	  InfiniteVector<double, Index> transformcoeff;
	  basis.reconstruct_t(wcoeff, level, transformcoeff);
	  
	  cout << "* generator: " << index
	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == last_generator(&basis, level)) break;
	}
    }
#endif
}

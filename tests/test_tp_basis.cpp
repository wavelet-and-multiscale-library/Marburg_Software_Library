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

#if 1
  cout << "- testing iterator functionality:" << endl;
  for (Index index(first_generator(&basis, basis.j0()));; ++index) {
    cout << index << endl;
    if (index == last_wavelet(&basis, basis.j0()+1)) break;
  }
#endif

#if 0
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(first_generator(&basis, level));
      for (;; ++index)
	{
	  cout << "* generator: " << index << endl;

// 	  InfiniteVector<double, Index> origcoeff;
// 	  origcoeff[index] = 1.0;
	  
// 	  InfiniteVector<double, Index> wcoeff;
// 	  basis.decompose(origcoeff, basis.j0(), wcoeff);
	  
// 	  InfiniteVector<double, Index> transformcoeff;
// 	  basis.reconstruct(wcoeff, level, transformcoeff);
	  
// 	  cout << "* generator: " << index
// 	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == last_generator(&basis, level)) break;
	}
    }
#endif
  
}

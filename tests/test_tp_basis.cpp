#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>

#include <interval/spline_basis.h>
#include <generic/tp_basis.h>
#include <Rd/haar_mask.h>
#include <Rd/r_basis.h>
#include <interval/periodic.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing tensor product wavelet bases..." << endl;

#if 1
  typedef SplineBasis<2,2,P_construction,0,0,0,0,SplineBasisData_j0<2,2,P_construction,0,0,0,0>::j0> Basis1;
  typedef Basis1::Index Index1;

  typedef SplineBasis<2,2,P_construction,1,0,0,0,SplineBasisData_j0<2,2,P_construction,1,0,0,0>::j0> Basis2;
  typedef Basis2::Index Index2;

  typedef TensorProductBasis<Basis1,Basis2> Basis;
  typedef Basis::Index Index;
  Basis basis;
  
  cout << "* a tensor product of two PBasis<2,2> bases:" << endl;
#else
  typedef PeriodicBasis<RBasis<HaarMask> > Basis1;
  typedef Basis1::Index Index1;

  typedef PeriodicBasis<RBasis<HaarMask> > Basis2;
  typedef Basis2::Index Index2;

  typedef TensorProductBasis<Basis1,Basis2> Basis;
  typedef Basis::Index Index;
  Basis basis;

  cout << "* a tensor product of 2 Haar bases:" << endl;
#endif

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- first generator on the coarsest level: " << basis.first_generator(basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << basis.last_generator(basis.j0()) << endl;
  cout << "- first wavelet on the coarsest level: " << basis.first_wavelet(basis.j0()) << endl;
  cout << "- last wavelet on the coarsest level: " << basis.last_wavelet(basis.j0()) << endl;

#if 1
  cout << "- testing iterator functionality:" << endl;
  for (Index index(basis.first_generator(basis.j0()));; ++index) {
    cout << index << endl;
    if (index == basis.last_wavelet(basis.j0()+1)) break;
  }
#endif

#if 1
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(basis.first_generator(level));
      for (;; ++index)
	{
 	  InfiniteVector<double, Index> origcoeff;
 	  origcoeff[index] = 1.0;
	  
// 	  cout << "* original coeffs:" << endl << origcoeff;

 	  InfiniteVector<double, Index> wcoeff;
 	  basis.decompose(origcoeff, basis.j0(), wcoeff);
	  
// 	  cout << "* after decompose():" << endl << wcoeff;

 	  InfiniteVector<double, Index> transformcoeff;
 	  basis.reconstruct(wcoeff, level, transformcoeff);
	  
// 	  cout << "* after reconstruct():" << endl << transformcoeff;

 	  cout << "* generator: " << index
 	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == basis.last_generator(level)) break;
	}
    }
#endif
  
#if 0
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose_t() and reconstruct_t() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(basis.first_generator(level));
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
	  
	  if (index == basis.last_generator(level)) break;
	}
    }
#endif
}

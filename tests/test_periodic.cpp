#include <iostream>
#include <Rd/daubechies_mask.h>
#include <Rd/r_basis.h>
#include <Rd/cdf_basis.h>
#include <interval/periodic.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing periodic wavelet bases ..." << endl;

  typedef PeriodicBasis<RBasis<DaubechiesMask<3> > > Basis;
  typedef Basis::Index Index;
  Basis basis;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;

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

	  InfiniteVector<double, Index> wcoeff;
	  basis.decompose(origcoeff, basis.j0(), wcoeff);

	  InfiniteVector<double, Index> transformcoeff;
	  basis.reconstruct(wcoeff, level, transformcoeff);

	  cout << "* generator: " << index
	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == last_generator(&basis, level)) break;
	}
    }
#endif

#if 1
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


  // ---------------------------------------------------------------------------

  cout << "* another periodic basis:" << endl;

  typedef PeriodicBasis<CDFBasis<2, 2> > Basis2;
  typedef Basis2::Index Index2;
  Basis2 basis2;

  cout << "- j0=" << basis2.j0() << endl;
  cout << "- the default wavelet index: " << Index2() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis2, basis2.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis2, basis2.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis2, basis2.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis2, basis2.j0()) << endl;

#if 1
  for (int level = basis2.j0()+1; level <= basis2.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index2 index(first_generator(&basis2, level));
      for (;; ++index)
	{
	  InfiniteVector<double, Index2> origcoeff;
	  origcoeff[index] = 1.0;

	  InfiniteVector<double, Index2> wcoeff;
	  basis2.decompose(origcoeff, basis2.j0(), wcoeff);

	  InfiniteVector<double, Index2> transformcoeff;
	  basis2.reconstruct(wcoeff, level, transformcoeff);

	  cout << "* generator: " << index
	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == last_generator(&basis2, level)) break;
	}
    }
#endif

#if 1
  for (int level = basis2.j0()+1; level <= basis2.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index2 index(first_generator(&basis2, level));
      for (;; ++index)
	{
	  InfiniteVector<double, Index2> origcoeff;
	  origcoeff[index] = 1.0;
	  
	  InfiniteVector<double, Index2> wcoeff;
 	  basis2.decompose(origcoeff, basis2.j0(), wcoeff);
	  
	  InfiniteVector<double, Index2> transformcoeff;
	  basis2.reconstruct(wcoeff, level, transformcoeff);

	  double error = linfty_norm(origcoeff-transformcoeff);
	  
	  cout << "* generator: " << index
	       << ", max. error: " << error << endl;

	  if (index == last_generator(&basis2, level)) break;
	}
    }
#endif

  return 0;
}

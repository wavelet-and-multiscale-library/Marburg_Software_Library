#include <iostream>
#include <fstream>
#include <Rd/r_basis.h>
#include <Rd/cdf_basis.h>
#include <interval/periodic.h>
#include <interval/i_indexplot.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing periodic wavelet bases ..." << endl;

  cout << "* a periodic Haar basis:" << endl;

  typedef PeriodicBasis<CDFBasis<1, 1> > Basis;
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
      Index index(first_generator<Basis>(level));
      for (;; ++index)
 	{
 	  InfiniteVector<double, Index> origcoeff;
 	  origcoeff[index] = 1.0;

 	  InfiniteVector<double, Index> wcoeff;
 	  basis.decompose(origcoeff, Basis::j0(), wcoeff);

 	  InfiniteVector<double, Index> transformcoeff;
 	  basis.reconstruct(wcoeff, level, transformcoeff);

 	  cout << "* generator: " << index
 	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
 	  if (index == last_generator<Basis>(level)) break;
 	}
    }
#endif

#if 1
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose_t() and reconstruct_t() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(first_generator<Basis>(level));
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
	  
	  if (index == last_generator<Basis>(level)) break;
	}
    }
#endif


#if 1
  // ---------------------------------------------------------------------------

  cout << "* a periodic CDF basis:" << endl;

  typedef PeriodicBasis<CDFBasis<2, 2> > Basis2;
  typedef Basis2::Index Index2;
  Basis2 basis2;

  cout << "- j0=" << basis2.j0() << endl;
  cout << "- the default wavelet index: " << Index2() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis2, basis2.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis2, basis2.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis2, basis2.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis2, basis2.j0()) << endl;

  for (int level = basis2.j0()+1; level <= basis2.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index2 index(first_generator<Basis2>(level));
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
	  
	  if (index == last_generator<Basis2>(level)) break;
	}
    }

  for (int level = basis2.j0()+1; level <= basis2.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index2 index(first_generator<Basis2>(level));
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

	  if (index == last_generator<Basis2>(level)) break;
	}
    }
#endif

#if 1
  cout << "- create some test index set..." << endl;
  InfiniteVector<double, Index2> gcoeffs, coeffs;
  gcoeffs[++(++first_generator<Basis2>(basis2.j0()+3))] = 1.0;
  cout << "* original coefficient set:" << endl
       << gcoeffs << endl;

  basis2.decompose(gcoeffs, basis2.j0(), coeffs);
  cout << "* coefficient set in multiscale representation:" << endl
       << coeffs << endl;
  
  cout << "* plotting the coefficient set..." << endl;
  std::ofstream fs("coefficient_plot.m");
  plot_indices(&basis2, coeffs, basis2.j0()+2, fs);
  fs.close();
  cout << "  ...done!" << endl;
#endif

  return 0;
}

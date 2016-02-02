#include <iostream>
#include <fstream>
#include <utils/array1d.h>
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

#if 0
  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(basis.first_generator(level));
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

#if 1
  {
    const int j0 = basis.j0();
    Vector<double> x(basis.Deltasize(j0+1));
    x[0] = 1;
    cout << "* a vector x=" << x << endl;
    Vector<double> y(basis.Deltasize(j0+1));
    basis.apply_Mj(j0, x, y);
    cout << "* applying Mj=(Mj0 Mj1) to x yields y=Mj*x=" << y << endl;
    basis.apply_Gj(j0, y, x);
    cout << "* applying Gj=(Mj0T Mj1T)^T to y yields x=Gj*y=" << x << endl;
    basis.apply_Gj_transposed(j0, x, y);
    cout << "* applying Gj^T to x yields y=Gj^T*x=" << y << endl;
    basis.apply_Mj_transposed(j0, x, y);
    cout << "* applying Mj^T to y yields x=Mj^T*y=" << x << endl;

    x.scale(0); y.scale(0);
    x[0] = 1;
    basis.apply_Tj(j0, x, y);
    cout << "* applying T_{j0} to x yields y=" << y << endl;
    x.resize(basis.Deltasize(j0+2));
    x[3] = 1;
    cout << "* x on the next level: " << x << endl;
    y.resize(basis.Deltasize(j0+2));
    basis.apply_Tj(j0+1, x, y);
    cout << "* applying T_{j0+1} to x yields y=" << y << endl;
    basis.apply_Tjinv(j0+1, y, x);
    cout << "* applying T_{j0+1}^{-1} to y yields x=" << x << endl;
    x.resize(basis.Deltasize(j0+3));
    x[1] = 1;
    cout << "* x on the next plus 1 level: " << x << endl;
    y.resize(basis.Deltasize(j0+3));
    basis.apply_Tj(j0+2, x, y);
    cout << "* applying T_{j0+2} to x yields y=" << y << endl;
    basis.apply_Tjinv(j0+2, y, x);
    cout << "* applying T_{j0+2}^{-1} to y yields x=" << x << endl;
  }
#endif


  // ---------------------------------------------------------------------------

  cout << "* a periodic CDF basis:" << endl;

  typedef PeriodicBasis<CDFBasis<2,2> > Basis2;
  typedef Basis2::Index Index2;
  Basis2 basis2;

  cout << "- j0=" << basis2.j0() << endl;
  cout << "- the default wavelet index: " << Index2() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis2, basis2.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis2, basis2.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis2, basis2.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis2, basis2.j0()) << endl;

#if 0
  for (int level = basis2.j0()+1; level <= basis2.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index2 index(basis2.first_generator(level));
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
	  
	  if (index == basis2.last_generator(level)) break;
	}
    }

  for (int level = basis2.j0()+1; level <= basis2.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index2 index(basis2.first_generator(level));
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

	  if (index == basis2.last_generator(level)) break;
	}
    }
#endif

#if 0
  cout << "- create some test index set..." << endl;
  InfiniteVector<double, Index2> gcoeffs, coeffs;
  gcoeffs[++(++basis2.first_generator(basis2.j0()+3))] = 1.0;
  cout << "* original coefficient set:" << endl
       << gcoeffs << endl;

  basis2.decompose(gcoeffs, basis2.j0(), coeffs);
  cout << "* coefficient set in multiscale representation:" << endl
       << coeffs << endl;
  
  cout << "* plotting the coefficient set..." << endl;
  std::ofstream fs("coefficient_plot.m");
  plot_indices2(&basis2, coeffs, basis2.j0()+2, fs);
  fs.close();
  cout << "  ...done!" << endl;
#endif

#if 0
  cout << "* point evaluation of periodic CDF functions:" << endl;
  int N = 32;
  Array1D<double> points(N+1), values(N+1), dervalues(N+1);
  double h = 1.0/N;
  for (int i = 0; i <= N; i++) points[i] = i*h;
  const int level = basis2.j0();
  for (Index2 lambda(basis2.first_generator(level));; ++lambda) {
    cout << lambda << endl;
    basis2.evaluate(lambda, points, values, dervalues);
    cout << "points: " << points << endl;
    cout << "values: " << values << endl;
    cout << "values of first derivative: " << dervalues << endl;
    if (lambda == basis2.last_wavelet(level)) break;
  }
#endif

#if 1
  {
    const int j0 = basis2.j0();
    Vector<double> x(basis2.Deltasize(j0+1));
    x[0] = 1;
    cout << "* a vector x=" << x << endl;
    Vector<double> y(basis2.Deltasize(j0+1));
    basis2.apply_Mj(j0, x, y);
    cout << "* applying Mj=(Mj0 Mj1) to x yields y=Mj*x=" << y << endl;
    basis2.apply_Gj(j0, y, x);
    cout << "* applying Gj=(Mj0T Mj1T)^T to y yields x=Gj*y=" << x << endl;
    basis2.apply_Gj_transposed(j0, x, y);
    cout << "* applying Gj^T to x yields y=Gj^T*x=" << y << endl;
    basis2.apply_Mj_transposed(j0, x, y);
    cout << "* applying Mj^T to y yields x=Mj^T*y=" << x << endl;

    x.scale(0); y.scale(0);
    x[0] = 1;
    basis2.apply_Tj(j0, x, y);
    cout << "* applying T_{j0} to x yields y=" << y << endl;
    x.resize(basis2.Deltasize(j0+2));
    x[3] = 1;
    cout << "* x on the next level: " << x << endl;
    y.resize(basis2.Deltasize(j0+2));
    basis2.apply_Tj(j0+1, x, y);
    cout << "* applying T_{j0+1} to x yields y=" << y << endl;
    basis2.apply_Tjinv(j0+1, y, x);
    cout << "* applying T_{j0+1}^{-1} to y yields x=" << x << endl;
    x.resize(basis2.Deltasize(j0+3));
    x[1] = 1;
    cout << "* x on the next plus 1 level: " << x << endl;
    y.resize(basis2.Deltasize(j0+3));
    basis2.apply_Tj(j0+2, x, y);
    cout << "* applying T_{j0+2} to x yields y=" << y << endl;
    basis2.apply_Tjinv(j0+2, y, x);
    cout << "* applying T_{j0+2}^{-1} to y yields x=" << x << endl;
  }
#endif

  return 0;
}

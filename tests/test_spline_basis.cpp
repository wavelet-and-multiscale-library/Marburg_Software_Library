#include <iostream>
#include <fstream>
#include <algebra/vector.h>
#include <interval/spline_basis_data.h>
#include <interval/spline_basis.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing setup of SplineBasisData objects..." << endl;

#if 0
  // PBasis, no b.c.'s
  SplineBasisData<2,2,P_construction,0,0,0,0> sd22nobc(SplineBasisData_j0<2,2,P_construction,0,0,0,0>::j0);
  sd22nobc.check();
  
  // PBasis, complementary b.c.'s
  SplineBasisData<2,2,P_construction,1,1,0,0> sd22(SplineBasisData_j0<2,2,P_construction,1,1,0,0>::j0); 
  sd22.check();
  
  // PBasis, complementary b.c.'s at x=0
  SplineBasisData<2,2,P_construction,1,0,0,0> sd22bcleft(SplineBasisData_j0<2,2,P_construction,1,0,0,0>::j0);
  sd22bcleft.check();
  
  // PBasis, complementary b.c.'s at x=1
  SplineBasisData<2,2,P_construction,0,1,0,0> sd22bcright(SplineBasisData_j0<2,2,P_construction,0,1,0,0>::j0);
  sd22bcright.check();
  
  // PBasis, complementary b.c.'s at x=0
  SplineBasisData<3,3,P_construction,1,0,0,0> sd33bcleft(SplineBasisData_j0<3,3,P_construction,1,0,0,0>::j0);
  sd33bcleft.check();

  // PBasis, complementary b.c.'s at x=1
  SplineBasisData<3,3,P_construction,0,1,0,0> sd33bcright(SplineBasisData_j0<3,3,P_construction,0,1,0,0>::j0); 
  sd33bcright.check();
  
  // PBasis, complementary b.c.'s
  SplineBasisData<3,3,P_construction,1,1,0,0> sd33(SplineBasisData_j0<3,3,P_construction,1,1,0,0>::j0);
  sd33.check();
#endif
 
#if 0
  // DSBasis, no b.c.'s
  SplineBasisData<2,2,DS_construction_bio5,0,0,0,0> sb22nobc(SplineBasisData_j0<2,2,DS_construction_bio5,0,0,0,0>::j0);
  sb22nobc.check();
  
  // DSBasis, no b.c.'s
  SplineBasisData<2,2,DS_construction_bio5e,0,0,0,0> sb22nobc_energy(SplineBasisData_j0<2,2,DS_construction_bio5e,0,0,0,0>::j0);
  sb22nobc_energy.check();
  
  // DSBasis, no b.c.'s
  SplineBasisData<3,3,DS_construction_bio5,0,0,0,0> sb33nobc(SplineBasisData_j0<3,3,DS_construction_bio5,0,0,0,0>::j0);
  sb33nobc.check();
  
  // DSBasis, no b.c.'s
  SplineBasisData<3,3,DS_construction_bio5e,0,0,0,0> sb33nobc_energy(SplineBasisData_j0<3,3,DS_construction_bio5e,0,0,0,0>::j0); 
  sb33nobc_energy.check();
#endif

#if 1
  cout << "Testing SplineBasis..." << endl;
  // PBasis, complementary b.c.'s
  typedef SplineBasis<3,3,P_construction,1,1,0,0,SplineBasisData_j0<3,3,P_construction,1,1,0,0>::j0> SBasis;
  //   typedef SplineBasis<3,3,DS_construction_bio5e,0,0,0,0> SBasis; // DSBasis, no b.c.'s
  SBasis basis;
  typedef SBasis::Index Index;

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
  x[4] = 1;
  cout << "* x on the next level: " << x << endl;
  y.resize(basis.Deltasize(j0+2));
  basis.apply_Tj(j0+1, x, y);
  cout << "* applying T_{j0+1} to x yields y=" << y << endl;
  basis.apply_Tjinv(j0+1, y, x);
  cout << "* applying T_{j0+1}^{-1} to y yields x=" << x << endl;
  x.resize(basis.Deltasize(j0+3));
  x[4] = 1;
  cout << "* x on the next plus 1 level: " << x << endl;
  y.resize(basis.Deltasize(j0+3));
  basis.apply_Tj(j0+2, x, y);
  cout << "* applying T_{j0+2} to x yields y=" << y << endl;
  basis.apply_Tjinv(j0+2, y, x);
  cout << "* applying T_{j0+2}^{-1} to y yields x=" << x << endl;
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
	  
  	  cout << "* original coeffs:" << endl << origcoeff;

	  InfiniteVector<double, Index> wcoeff;
	  basis.decompose(origcoeff, basis.j0(), wcoeff);

  	  cout << "* after decompose():" << endl << wcoeff;
	  
	  InfiniteVector<double, Index> transformcoeff;
	  basis.reconstruct(wcoeff, level, transformcoeff);

	  cout << "* after reconstruct():" << endl << transformcoeff;
	  
	  cout << "* generator: " << index
	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == basis.last_generator(level)) break;
	}
    }
#endif

#if 0  
  cout << "* point evaluation of spline wavelets:" << endl;
  typedef SBasis::Index Index;
  int N = 32;
  Array1D<double> points(N+1), values(N+1), dervalues(N+1);
  double h = 1.0/N;
  for (int i = 0; i <= N; i++) points[i] = i*h;
  const int level = basis.j0();
  for (Index lambda(basis.first_generator(level));; ++lambda) {
    cout << lambda << endl;
    basis.evaluate(lambda, points, values, dervalues);
    cout << "values: " << values << endl;
    cout << "values of first derivative: " << dervalues << endl;
//     if (lambda == basis.last_generator(level)) break;
    if (lambda == basis.last_wavelet(level)) break;
  }
#endif

  return 0;
}

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

#if 0
  typedef SplineBasis<2,2,P_construction,0,0,0,0,SplineBasisData_j0<2,2,P_construction,0,0,0,0>::j0> Basis0;
  typedef Basis0::Index Index0;

  typedef SplineBasis<2,2,P_construction,1,0,0,0,SplineBasisData_j0<2,2,P_construction,1,0,0,0>::j0> Basis1;
  typedef Basis1::Index Index1;

  typedef TensorProductBasis<Basis0,Basis1> Basis;
  typedef Basis::Index Index;
  Basis basis;
  
  cout << "* a tensor product of two PBasis<2,2> bases:" << endl;
#else
  typedef PeriodicBasis<RBasis<HaarMask> > Basis0;
  typedef Basis0::Index Index0;

  typedef PeriodicBasis<RBasis<HaarMask> > Basis1;
  typedef Basis1::Index Index1;

  typedef TensorProductBasis<Basis0,Basis1> Basis;
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

#if 1
  {
    const int j0 = basis.j0();
    Vector<double> x(basis.Deltasize(j0+1));
    x[3] = 1;
    cout << "* a vector x=" << x << endl;
    Vector<double> y(basis.Deltasize(j0+1));
    basis.apply_Mj(j0, x, y);
    cout << "* applying Mj=(Mj0 Mj1) to x yields y=Mj*x=" << y << endl;
    basis.apply_Gj(j0, y, x);
    cout << "* applying Gj=(Mj0T Mj1T)^T to y yields x=Gj*y=" << x << endl;
//     basis.apply_Gj_transposed(j0, x, y);
//     cout << "* applying Gj^T to x yields y=Gj^T*x=" << y << endl;
//     basis.apply_Mj_transposed(j0, x, y);
//     cout << "* applying Mj^T to y yields x=Mj^T*y=" << x << endl;

//     x.scale(0); y.scale(0);
//     x[0] = 1;
//     basis.apply_Tj(j0, x, y);
//     cout << "* applying T_{j0} to x yields y=" << y << endl;
//     x.resize(basis.Deltasize(j0+2));
//     x[3] = 1;
//     cout << "* x on the next level: " << x << endl;
//     y.resize(basis.Deltasize(j0+2));
//     basis.apply_Tj(j0+1, x, y);
//     cout << "* applying T_{j0+1} to x yields y=" << y << endl;
//     basis.apply_Tjinv(j0+1, y, x);
//     cout << "* applying T_{j0+1}^{-1} to y yields x=" << x << endl;
//     x.resize(basis.Deltasize(j0+3));
//     x[1] = 1;
//     cout << "* x on the next plus 1 level: " << x << endl;
//     y.resize(basis.Deltasize(j0+3));
//     basis.apply_Tj(j0+2, x, y);
//     cout << "* applying T_{j0+2} to x yields y=" << y << endl;
//     basis.apply_Tjinv(j0+2, y, x);
//     cout << "* applying T_{j0+2}^{-1} to y yields x=" << x << endl;
  }
#endif

}

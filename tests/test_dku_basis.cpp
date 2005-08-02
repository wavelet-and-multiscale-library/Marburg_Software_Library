#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>

#include <interval/dku_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DKU bases..." << endl;

  const int d = 2;
  const int dT = 4;
  typedef DKUBasis<d, dT> Basis;
  Basis basis; // Bernstein SVD
//   Basis basis(none);
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;

  cout << "- the default wavelet index: " << Basis::Index(&basis) << endl;

  cout << "- leftmost generator on the coarsest level: " << basis.firstGenerator(basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << basis.lastGenerator(basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << basis.firstWavelet(basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << basis.lastWavelet(basis.j0()) << endl;

  typedef Basis::Index Index;

#if 0
  cout << "- iterating from first generator on coarsest level to last wavelet on next level:" << endl;
  Index index(basis.firstGenerator(basis.j0()));
  for (;; ++index) {
    cout << index << endl;
    if (index == basis.lastWavelet(basis.j0()+1)) break;
  }
#endif

  cout << "- index of leftmost right boundary generator on scale j0: "
       << basis.DeltaRmin(basis.j0()) << endl;

  cout << "- checking biorthogonality of Mj0, Mj0T for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++)
    {
      SparseMatrix<double> mj0_t, mj0T;
      basis.assemble_Mj0_t(level, mj0_t);
      basis.assemble_Mj0T(level, mj0T);

      SparseMatrix<double> T = mj0_t * mj0T;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj0^T*Mj0T-I||_infty: " << row_sum_norm(T) << endl;

      SparseMatrix<double> mj0, mj0T_t;
      basis.assemble_Mj0(level, mj0);
      basis.assemble_Mj0T_t(level, mj0T_t);

      T = mj0T_t * mj0;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj0T^T*Mj0-I||_infty: " << row_sum_norm(T) << endl;
    }

  cout << "- checking biorthogonality of Mj<->Gj and MjT<->GjT for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+1; level++)
    {
      SparseMatrix<double> mj0, mj1;
      basis.assemble_Mj0(level, mj0);
      basis.assemble_Mj1(level, mj1);
      SparseMatrix<double> mj(mj0.row_dimension(), mj0.row_dimension());
      mj.set_block(0, 0, mj0);
      mj.set_block(0, mj0.column_dimension(), mj1);

      SparseMatrix<double> mj0T_t, mj1T_t;
      basis.assemble_Mj0T_t(level, mj0T_t);
      basis.assemble_Mj1T_t(level, mj1T_t);
      SparseMatrix<double> gj(mj0T_t.column_dimension(), mj0T_t.column_dimension());
      gj.set_block(0, 0, mj0T_t);
      gj.set_block(mj0T_t.row_dimension(), 0, mj1T_t);

      SparseMatrix<double> T = mj * gj;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj*Gj-I||_infty: " << row_sum_norm(T) << endl;

      T = gj * mj;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Gj*Mj-I||_infty: " << row_sum_norm(T) << endl;
      
      SparseMatrix<double> mj0T, mj1T;
      basis.assemble_Mj0T(level, mj0T);
      basis.assemble_Mj1T(level, mj1T);
      SparseMatrix<double> mjt(mj.row_dimension(), mj.row_dimension());
      mjt.set_block(0, 0, mj0T);
      mjt.set_block(0, mj0T.column_dimension(), mj1T);

      SparseMatrix<double> mj0_t, mj1_t;
      basis.assemble_Mj0_t(level, mj0_t);
      basis.assemble_Mj1_t(level, mj1_t);
      SparseMatrix<double> gjt(mj.row_dimension(), mj.row_dimension());
      gjt.set_block(0, 0, mj0_t);
      gjt.set_block(mj0_t.row_dimension(), 0, mj1_t);

      T = mjt * gjt;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||MjT*GjT-I||_infty: " << row_sum_norm(T) << endl;

      T = gjt * mjt;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||GjT*MjT-I||_infty: " << row_sum_norm(T) << endl;
    }

  for (int level = basis.j0()+1; level <= basis.j0()+2; level++)
    {
      cout << "- checking decompose() and reconstruct() for some/all generators on the level "
	   << level << ":" << endl;
      Index index(basis.firstGenerator(level));
//      for (int i = 1; i <= 4; ++i, ++index)
      for (;; ++index)
	{
	  InfiniteVector<double, Index> origcoeff;
	  origcoeff[index] = 1.0;
	  
// 	  cout << "original index set:" << endl << origcoeff;
	  
	  InfiniteVector<double, Index> wcoeff;
	  basis.decompose(origcoeff, basis.j0(), wcoeff);
	  
// 	  cout << "wavelet coefficients:" << endl << wcoeff;
	  
	  InfiniteVector<double, Index> transformcoeff;
	  basis.reconstruct(wcoeff, level, transformcoeff);
	  
// 	  cout << "after decompose()+reconstruct():" << endl << origcoeff;
	  
	  cout << "* generator: " << index
	       << ", max. error: " << linfty_norm(origcoeff-transformcoeff) << endl;
	  
	  if (index == basis.lastGenerator(level)) break;
	}
    }

#if 0
//   coeff[basis.firstGenerator(basis.j0())] = 0.0;
//   coeff[++basis.firstGenerator(basis.j0())] = 0.0;
  Index index(basis.firstGenerator(basis.j0()+1));
  for (int i = 0; i <= 15; i++, ++index);
  coeff[index] = 1.0;
//   for (int i = 0; i <= 8; i++)
//     coeff[++index] = 1.0;
  cout << "  * a small but nontrivial coefficient set:" << endl << coeff;
  cout << "  * result of DECOMPOSE:" << endl;
  InfiniteVector<double,Index> wcoeff;
  basis.decompose(coeff, basis.j0(), wcoeff);
  cout << wcoeff;

  cout << "  * RECONSTRUCT that:" << endl;
  InfiniteVector<double,Index> rcoeff;
  basis.reconstruct(wcoeff,basis.j0()+1,rcoeff);
  cout << rcoeff;
#endif

#if 0
  cout << "- evaluating some primal generators:" << endl;
  Index lambda(basis.firstGenerator(basis.j0()));
  for (;; ++lambda) {
    cout << lambda << endl;
    basis.evaluate(lambda, true, 5).matlab_output(cout);
    if (lambda == basis.lastGenerator(basis.j0())) break;
  }
#endif

#if 0
  cout << "- evaluating some dual generators:" << endl;
  lambda = basis.firstGenerator(basis.j0());
  for (;; ++lambda) {
    basis.evaluate(lambda, false, 6).matlab_output(cout);
    if (lambda == basis.lastGenerator(basis.j0())) break;
  }
#endif

  cout << "* another basis:" << endl;

  const int d2 = 3;
  const int dT2 = 5;
  DKUBasis<d2, dT2> basis2;

  cout << "- the (" << d2 << "," << dT2 << ") basis has j0=" << basis2.j0() << endl;
  cout << "- the default wavelet index: " << DKUBasis<d2, dT2>::Index(&basis2) << endl;
  cout << "- leftmost generator on the coarsest level: " << basis2.firstGenerator(basis2.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << basis2.lastGenerator(basis2.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << basis2.firstWavelet(basis2.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << basis2.lastWavelet(basis2.j0()) << endl;

  cout << "* yet another basis:" << endl;

//   const int d3 = 3;
//   const int dT3 = 7;
//   DKUBasis<d3, dT3> basis3;

//   cout << "- the (" << d3 << "," << dT3 << ") basis has j0=" << basis3.j0() << endl;
//   cout << "- the default wavelet index: " << DKUBasis<d3, dT3>::Index(&basis3) << endl;
//   cout << "- leftmost generator on the coarsest level: " << basis3.firstGenerator(basis3.j0()) << endl;
//   cout << "- rightmost generator on the coarsest level: " << basis3.lastGenerator(basis3.j0()) << endl;
//   cout << "- leftmost wavelet on the coarsest level: " << basis3.firstWavelet(basis3.j0()) << endl;
//   cout << "- rightmost wavelet on the coarsest level: " << basis3.lastWavelet(basis3.j0()) << endl;

  return 0;
}

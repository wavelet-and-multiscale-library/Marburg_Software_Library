#include <iostream>

#include <algebra/infinite_vector.h>
#include <algebra/sparse_matrix.h>
#include <utils/array1d.h>

#include <Rd/cdf_utils.h>
#include <interval/ds_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DS bases..." << endl;

  const int d = 2;
  const int dT = 2;

  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;

  Basis basis; // Z={0,1}
//   Basis basis(1, 0, 0, 1); // Z={0}
//   Basis basis(0, 1, 1, 0); // Z={1}
//   Basis basis(0, 0, 1, 1); // Z={}
//   Basis basis(0, 0, 0, 0); // should work, DKU basis without b.c.'s at all
  
  cout << "- d=" << d << ", dT=" << dT << endl;
  cout << "- ell1=" << ell1<d>() << ", ell2=" << ell2<d>()
       << ", ell1T=" << ell1T<d,dT>() << ", ell2T=" << ell2T<d,dT>() << endl;
  cout << "- ellT_l=" << basis.ellT_l() << ", ellT_r=" << basis.ellT_r()
       << ", ell_l=" << basis.ell_l() << ", ell_r=" << basis.ell_r() << endl;

  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;

#if 1
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
#endif


  cout << "* another basis:" << endl;

  const int d2 = 3;
  const int dT2 = 5;
  
  typedef DSBasis<d2, dT2> Basis2;
  typedef Basis2::Index Index2;

//   Basis2 basis2; // Z={0,1}
//   Basis2 basis2(1, 0, 0, 1); // Z={0}
  Basis2 basis2(0, 1, 1, 0); // Z={1}
//   Basis2 basis2(0, 0, 1, 1); // Z={}
//   Basis2 basis2(0, 0, 0, 0); // should work, DKU basis without b.c.'s at all

  cout << "- d=" << d2 << ", dT=" << dT2 << endl;
  cout << "- ell1=" << ell1<d2>() << ", ell2=" << ell2<d2>()
       << ", ell1T=" << ell1T<d2,dT2>() << ", ell2T=" << ell2T<d2,dT2>() << endl;
  cout << "- ellT_l=" << basis2.ellT_l() << ", ellT_r=" << basis2.ellT_r()
       << ", ell_l=" << basis2.ell_l() << ", ell_r=" << basis2.ell_r() << endl;

  cout << "- the (" << d2 << "," << dT2 << ") basis has j0=" << basis2.j0() << endl;
  cout << "- the default wavelet index: " << Index2() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis2, basis2.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis2, basis2.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis2, basis2.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis2, basis2.j0()) << endl;

#if 1
  cout << "- checking biorthogonality of Mj0, Mj0T for different levels:" << endl;
  for (int level = basis2.j0(); level <= basis2.j0()+2; level++)
    {
      SparseMatrix<double> mj0_t, mj0T;
      basis2.assemble_Mj0_t(level, mj0_t);
      basis2.assemble_Mj0T(level, mj0T);

      SparseMatrix<double> T = mj0_t * mj0T;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj0^T*Mj0T-I||_infty: " << row_sum_norm(T) << endl;

      SparseMatrix<double> mj0, mj0T_t;
      basis2.assemble_Mj0(level, mj0);
      basis2.assemble_Mj0T_t(level, mj0T_t);

      T = mj0T_t * mj0;
      for (unsigned int i = 0; i < T.row_dimension(); i++)
	T.set_entry(i, i, T.get_entry(i, i) - 1.0);
      cout << "* j=" << level << ",  ||Mj0T^T*Mj0-I||_infty: " << row_sum_norm(T) << endl;
    }
#endif

}

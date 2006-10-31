#include <iostream>
#include <interval/p_basis.h>
#include <interval/ds_basis.h>
#include <numerics/bvp.h>
#include <cube/cube_basis.h>
#include <galerkin/cube_equation.h>

using WaveletTL::CubeEquation;
using WaveletTL::CubeBasis;
using MathTL::IdentityBVP;


using namespace std;
using namespace WaveletTL;
using namespace MathTL;

int main()
{
  cout << "Testing wavelet bases from [P] ..." << endl;

  const int d  = 3;
  const int dT = 3;


  typedef PBasis<d,dT> Basis;
  typedef Basis::Index Index;

  Basis basis; // no b.c.'s
//   Basis basis(1, 1); // 1st order complementary b.c.'s at x=0 and x=1
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1
//   Basis basis(3, 3);

  cout << "- d=" << d << ", dT=" << dT << endl;
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(&basis, basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(&basis, basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(&basis, basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(&basis, basis.j0()) << endl;


  //  abort();  

#if 0
  {
    SparseMatrix<double> mj0;
    const int j = basis.j0();
    basis.assemble_Mj0(j, mj0);
    cout << "* Mj0 for j=j0=" << j << ":" << endl;
    cout << mj0;
    std::ofstream result("Mj0.m");
    result << "M=";
    print_matrix(mj0, result);
    result << endl;
    result.close();
  }
#endif


#if 0
  
  typedef CubeBasis<Basis,1> CBasis;
  typedef CBasis::Index CIndex;

  FixedArray1D<int,2> bc;
  bc[0] = 1;
  bc[1] = 2;
  
  CBasis cbasis(bc);

  set<CIndex> Lambda;
  for (CIndex lambda(first_generator<Basis,1,CBasis>(&cbasis, cbasis.j0()));; ++lambda) {
    cout << lambda << endl;
    Lambda.insert(lambda);
    if (lambda == last_generator<Basis,1,CBasis>(&cbasis, cbasis.j0())) break;
  }

   Vector<double> value(1);
   value[0] = 1;
   ConstantFunction<1> const_fun(value);
   IdentityBVP<1> trivial_bvp(&const_fun);
   

   CubeEquation<Basis,1,CBasis> eq(&trivial_bvp, bc);


   cout << "setting up full right hand side..." << endl;
   Vector<double> rh;
   WaveletTL::setup_righthand_side(eq, Lambda, rh);
   cout << "...done setting up full right hand side" << endl;

   InfiniteVector<double, Index> coeff;
   Index index(first_generator(&basis, basis.j0()));
   for (int i = 0;; ++index, i++) {
     //cout << index << endl;
     coeff.set_coefficient(index, rh[i]);
     if (index == last_generator(&basis, basis.j0())) break;
   }
   
   //coeff.scale(&eq,-1);

   //cout << index << endl;
//    ++index;
//    ++index;
//    ++index;
//    ++index;
//    ++index;
//    coeff.set_coefficient(index, 1.);
   
   cout << "evaluating expansion in dual basis..." << endl;
   SampledMapping<1> res = evaluate(basis, coeff, false, 8);

   cout << "...done evaluating expansion in dual basis" << endl;

   
   std::ofstream ofs5("reproduced_function.m");
   res.matlab_output(ofs5);
   ofs5.close();

   //   abort();

#endif


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

#if 0
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
#endif

#if 0
  cout << "- checking access to single rows of the M_{j,i} matrices:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++)
    {
      InfiniteVector<double, Vector<double>::size_type> v, w;
      double maxerr = 0.0;
      SparseMatrix<double> mj0, mj0_t, mj1, mj1_t, mj0T, mj0T_t, mj1T, mj1T_t;
      basis.assemble_Mj0(level, mj0); mj0_t = transpose(mj0);
      basis.assemble_Mj1(level, mj1); mj1_t = transpose(mj1);
      basis.assemble_Mj0T(level, mj0T); mj0T_t = transpose(mj0T);
      basis.assemble_Mj1T(level, mj1T); mj1T_t = transpose(mj1T);
      for (size_t row = 0; row < mj0.row_dimension(); row++)
	{
	  mj0.get_row(row, v);
	  basis.Mj0_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj0_t.row_dimension(); row++)
	{
	  mj0_t.get_row(row, v);
	  basis.Mj0_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0_t: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1.row_dimension(); row++)
	{
	  mj1.get_row(row, v);
	  basis.Mj1_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj1: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1_t.row_dimension(); row++)
	{
	  mj1_t.get_row(row, v);
	  basis.Mj1_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj1_t: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj0T.row_dimension(); row++)
	{
	  mj0T.get_row(row, v);
	  basis.Mj0T_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0T: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj0T_t.row_dimension(); row++)
	{
	  mj0T_t.get_row(row, v);
	  basis.Mj0T_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj0T_t: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1.row_dimension(); row++)
	{
	  mj1T.get_row(row, v);
	  basis.Mj1T_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj1T: " << maxerr << endl;
      maxerr = 0.0;
      for (size_t row = 0; row < mj1T_t.row_dimension(); row++)
	{
	  mj1T_t.get_row(row, v);
	  basis.Mj1T_t_get_row(level, row, w);
	  maxerr = max(maxerr, linfty_norm(v-w));
	}
      cout << "* j=" << level << ", max. error in Mj1T_t: " << maxerr << endl;
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

#if 0
  cout << "- evaluating some primal generators:" << endl;
  for (Index lambda(first_generator(&basis, basis.j0()));; ++lambda) {
    cout << lambda << endl;
    evaluate(basis, lambda, true, 7).matlab_output(cout);
    if (lambda == last_generator(&basis, basis.j0())) break;
  }
  
  cout << "- evaluating some primal wavelets:" << endl;
  for (Index lambda = first_wavelet(&basis, basis.j0());; ++lambda) {
    cout << lambda << endl;
    evaluate(basis, lambda, true, 6).matlab_output(cout);
    if (lambda == last_wavelet(&basis, basis.j0())) break;
  }
#endif

  return 0;
}

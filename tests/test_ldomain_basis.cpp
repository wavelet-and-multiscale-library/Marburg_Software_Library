#include <iostream>
#include <map>

#include <algebra/infinite_vector.h>
#include <utils/map_tools.h>

#include <interval/ds_basis.h>
#include <interval/ds_evaluate.h>
#include <interval/spline_basis.h>

#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 2
#include <Ldomain/ldomain_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the L-shaped domain..." << endl;

  const int d  = 2;
  const int dT = 2;

#if 0
  typedef DSBasis<d,dT,BernsteinSVD> Basis1D;
  Basis1D basis1d;
#else
  typedef SplineBasis<d,dT,DS_construction> Basis1D;
  Basis1D basis1d("bio5-energy",0,0,0,0);
#endif
  
  typedef LDomainBasis<Basis1D> Basis;
  Basis basis(basis1d);

  typedef Basis::Index Index;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- first generator on the coarsest level: " << first_generator<Basis1D>(&basis, basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << last_generator<Basis1D>(&basis, basis.j0()) << endl;
  cout << "- first wavelet on the coarsest level: " << first_wavelet<Basis1D>(&basis, basis.j0()) << endl;
  cout << "- last wavelet on the coarsest level: " << last_wavelet<Basis1D>(&basis, basis.j0()) << endl;

#if 0
  // only for IBASIS != SplineBasis
  cout << "- checking setup of Mj0 for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++) {
    cout << "* j=" << level << endl;
    const BlockMatrix<double>& Mj0 = basis.get_Mj0(level); // should yield a cache miss
//     cout << "* j=" << level << ", Mj0=" << endl << Mj0 << endl;
    const BlockMatrix<double>& dummy = basis.get_Mj0(level); // should yield a cache hit
  } 

  cout << "- checking setup of Mj0T for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++) {
    cout << "* j=" << level << endl;
    const BlockMatrix<double>& Mj0T = basis.get_Mj0T(level); // should yield a cache miss
//     cout << "* j=" << level << ", Mj0T=" << endl << Mj0T << endl;
    const BlockMatrix<double>& dummy = basis.get_Mj0T(level); // should yield a cache hit
  } 
#endif

#if 0
  // only for IBASIS != SplineBasis
  cout << "- checking setup of Mj1c_1d for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++) {
    cout << "* j=" << level << endl;
    const SparseMatrix<double>& Mj1c_1d = basis.get_Mj1c_1d(level); // should yield a cache miss
//     cout << "* j=" << level << ", Mj1c_1d=" << endl << Mj1c_1d << endl;
    const SparseMatrix<double>& dummy = basis.get_Mj1c_1d(level); // should yield a cache hit
  } 
#endif

#if 0
  // only for IBASIS != SplineBasis
  cout << "- checking setup of Mj1c_01 for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++) {
    cout << "* j=" << level << endl;
    const BlockMatrix<double>& Mj1c_01 = basis.get_Mj1c_01(level); // should yield a cache miss
    cout << "* j=" << level << ", Mj1c_01=" << endl << Mj1c_01 << endl;
    const BlockMatrix<double>& dummy = basis.get_Mj1c_01(level); // should yield a cache hit
  } 
#endif

#if 0
  // only for IBASIS != SplineBasis
  cout << "- checking setup of Mj1c_10 for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++) {
    cout << "* j=" << level << endl;
    const BlockMatrix<double>& Mj1c_10 = basis.get_Mj1c_10(level); // should yield a cache miss
//     cout << "* j=" << level << ", Mj1c_10=" << endl << Mj1c_10 << endl;
    const BlockMatrix<double>& dummy = basis.get_Mj1c_10(level); // should yield a cache hit
  } 
#endif

#if 0
  // only for IBASIS != SplineBasis
  cout << "- checking setup of Mj1c_11 for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++) {
    cout << "* j=" << level << endl;
    const BlockMatrix<double>& Mj1c_11 = basis.get_Mj1c_11(level); // should yield a cache miss
//     cout << "* j=" << level << ", Mj1c_11=" << endl << Mj1c_11 << endl;
    const BlockMatrix<double>& dummy = basis.get_Mj1c_11(level); // should yield a cache hit
  } 
#endif

#if 0
  cout << "- checking index set sizes for different levels:" << endl;
  for (int level = basis.j0(); level <= basis.j0()+2; level++) {
    cout << "*   Deltasize(" << level << ")=" << basis.Deltasize(level) << endl;
    cout << "  Nabla01size(" << level << ")=" << basis.Nabla01size(level) << endl;
    cout << "  Nabla10size(" << level << ")=" << basis.Nabla10size(level) << endl;
    cout << "  Nabla11size(" << level << ")=" << basis.Nabla11size(level) << endl;
    cout << "  ------------------" << endl;
    cout << "  sum           ="
	 << basis.Deltasize(level)+basis.Nabla01size(level)+basis.Nabla10size(level)+basis.Nabla11size(level)
	 << endl;
  }
#endif

#if 1
 {
   cout << "- check apply_Mj0() and apply_Mj0T_transposed:" << endl;
   const int level = basis.j0();

   for (int i = 0; i < basis.Deltasize(level); i++) {
     std::map<size_t,double> x,y,z;
     x[i] = 1.0;
     basis.apply_Mj0(level,x,y);
     cout << "i=" << i << ": ";
//      cout << "Mj0*e_i=" << endl;
//      for (std::map<size_t,double>::const_iterator it(y.begin()); it != y.end(); ++it)
//        cout << it->first << ": " << it->second << endl;

     basis.apply_Mj0T_transposed(level,y,z);
//      cout << "Mj0T^T*Mj0*e_i=" << endl;
//      for (std::map<size_t,double>::const_iterator it(z.begin()); it != z.end(); ++it)
//        cout << it->first << ": " << it->second << endl;

     y.clear();
     add_maps(x,z,y,1.0,-1.0);

     cout << "||(I-Mj0T^T*Mj0)*e_i||_infty=";
     double error = 0;
     for (std::map<size_t,double>::const_iterator it(y.begin()); it != y.end(); ++it) {
//        if ((int)it->first < 3*(basis1d.Deltasize(basis.j0())-2)*(basis1d.Deltasize(basis.j0())-2))
//        if ((int)it->first >= 3*(basis1d.Deltasize(basis.j0())-2)*(basis1d.Deltasize(basis.j0())-2)
// 	   && (int)it->first < (3*(basis1d.Deltasize(basis.j0())-2)+1)*(basis1d.Deltasize(basis.j0())-2))
//        if ((int)it->first >= (3*(basis1d.Deltasize(basis.j0())-2)+1)*(basis1d.Deltasize(basis.j0())-2))
	 error = max(error, fabs(it->second));
     }
     cout << error << endl;
   }
 }
#endif

#if 0
  Index lambda(basis.first_generator(basis.j0()));
//   Index lambda(first_wavelet<Basis1D>(&basis, basis.j0()));
//   Index lambda(2, MultiIndex<int,2>(0,1), 0, MultiIndex<int,2>(3,0), &basis);

//   for (; !(lambda.p() == 1); ++lambda);
//   for (; !(lambda.p() == 2); ++lambda);
  for (; !(lambda.p() == 3); ++lambda);
//   for (; !(lambda.p() == 4); ++lambda);

//   for (; lambda.e()[0] != 0 || lambda.e()[1] != 1; ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 1); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 1 && lambda.k()[1] == 7); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 2); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 2 && lambda.k()[1] == 7); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 4); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 4 && lambda.k()[1] == 6); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 4 && lambda.k()[1] == 7); ++lambda);

//   for (; lambda.e()[0] != 1 || lambda.e()[1] != 0; ++lambda);

//   for (int i = 0; i < 155; i++, ++lambda); // one of the generators on patch 4
//   for (int i = 0; i < 330; i++, ++lambda); // one of the (0,1)-wavelets on patch 4
//   for (int i = 0; i < 334; i++, ++lambda);
  for (int i = 0; i < 2; i++, ++lambda);

  cout << "- evaluating a primal wavelet lambda=" << lambda << " ..." << endl;
  std::ofstream psistream("Ldomain_wavelet.m");
  matlab_output(psistream, basis.evaluate(lambda, 6));
//   octave_output(psistream, basis.evaluate(lambda, 6));
  psistream.close();
  cout << "  ...done, see file Ldomain_wavelet.m!" << endl;

  InfiniteVector<double, Index> gcoeffs;
  basis.reconstruct_1(lambda, lambda.j()+1, gcoeffs);
  cout << "- generator coefficients of lambda=" << lambda << " on a higher scale " << lambda.j()+1 << ":"
       << endl << gcoeffs;

  const int ecode = lambda.e()[0]+lambda.e()[1];
  if (ecode == 0) {
    std::map<size_t,double> gmap;
    gmap[lambda.number()] = 1.0;
//     for (std::map<size_t,double>::const_iterator it(gmap.begin()); it != gmap.end(); ++it)
//       cout << "gmap[" << it->first << "]=" << it->second << endl;
    
    std::map<size_t,double> gmap2;
    basis.apply_Mj0(lambda.j(), gmap, gmap2);
//     for (std::map<size_t,double>::const_iterator it(gmap2.begin()); it != gmap2.end(); ++it)
//       cout << "gmap2[" << it->first << "]=" << it->second << endl;

    InfiniteVector<double,Index> gcoeffs2;
    basis.map_to_vector(lambda.j()+1, gmap2, gcoeffs2);
//     cout << "- alternative computation of the generator coeffs yielded" << endl << gcoeffs2;
    cout << "- alternative computation of the generator coeffs yielded abs. error "
	 << linfty_norm(gcoeffs-gcoeffs2) << endl;
  }



  cout << "- evaluating this linear combination..." << endl;
  std::ofstream psi2stream("Ldomain_wavelet2.m");
  matlab_output(psi2stream, basis.evaluate(gcoeffs, 6));
  psi2stream.close();
  cout << "  ...done, see file Ldomain_wavelet2.m!" << endl;

#endif


  return 0;
}

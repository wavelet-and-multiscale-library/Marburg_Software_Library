#include <iostream>

#include <interval/p_basis.h>
#include <interval/ds_basis.h>
#include <interval/p_support.h>
//#include <cube/tbasis.h>
//#include <cube/tbasis_index.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
//#include <cube/cube_support.h>
#include <cmath>
#include <time.h> 

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of [P] generators and wavelets..." << endl;

  const int d = 3;
  const int dT = 3;
  const int dim = 2;

  typedef PBasis<d,dT> Basis1D;
  typedef Basis1D::Index Index1D;
//  typedef Basis1D::Support Support;
  typedef CubeBasis<Basis1D,dim> CBasis;
  //typedef TensorBasis<Basis1D,dim> TBasis;
  //typedef TBasis::Index Index;
//  typedef CBasis::Index Index;

  //! type index type
  typedef MultiIndex<int,dim> type_type;
  
  //! translation index type
  typedef MultiIndex<int,dim> translation_type;

//  typedef MultiIndex<int,dim> level_type;

//   Basis basis(0, 0); // no b.c.'s
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1

  FixedArray1D<int,2*dim> bound_1_2D;
  bound_1_2D[0] = 0;
  bound_1_2D[1] = 1;
  bound_1_2D[2] = 0;
  bound_1_2D[3] = 1;
  //bound_1_2D[4] = 0;
  //bound_1_2D[5] = 0;

  Basis1D basis1D(1,1);
  //CBasis basis(bound_1_2D); // complementary b.c.'s
  //TBasis basis(bound_1_2D); // complementary b.c.'s

  //Index lambda(first_generator<Basis1D,2,CBasis>(&basis, basis.j0()));
  //Index lambda(1, &basis);
  //Index1D lambda1D(0, &basis1D);
  //int j = 4;
  //int e = 0;
  //int k = 0;
  //level_type j;
//  int j;
  type_type e;
  translation_type k;
  e = type_type();
  //e[0]=1;
  //e[1]=1;
  //e[2]=1;

  clock_t tstart, tend;
//  double time;

  bool richtig = true;
  tstart = clock();
  Index1D lambda2(0, &basis1D);

  for( int i = 0 ; i<=70000 ; i++){
    //Index lambda(i, &basis);
    Index1D lambda(i, &basis1D);
    //lambda.get_Index_Parameter(&basis, i , j, e, k);
    //lambda1D.get_Index_Parameter(&basis1D, i , j, e, k);
    //cout << "Deltasize = " << basis.Deltasize(i) << "Nablasize = " << basis. Nablasize(i) << endl;
    //cout << "num =" << i  << "  lambda = " << lambda << " j = " << j << " e = " << e << " k = " << k  << endl;
    //cout << "num =" << i  << "  lambda1 = " << lambda << "  lambda2 = " << lambda2  << endl;
    //if(lambda.j() != j || lambda.e() != e || lambda.k() != k)
    if(lambda != lambda2)
      richtig = false;
    ++lambda2;
  }
  cout << " Alles richtig = " << richtig << endl;
  
  tend = clock();
  cout << " Alt braucht = " << (tend- tstart)/(double)CLOCKS_PER_SEC << endl;

  tstart = clock();
/*
  for( int i = 0 ; i<=60000 ; i++){
    //Index lambda(i, &basis);
    lambda.get_Index_Parameter(&basis, i , j, e, k);
  }
  tend = clock();
  cout << " Neu braucht = " << (tend- tstart)/(double)CLOCKS_PER_SEC << endl;
  */
#if 0
  cout << "- calculating some support intersections:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0());; ++lambda)
    {
      Support supp;
      support(basis, lambda, supp.k1, supp.k2);
      cout << "psi_lambda, lambda=" << lambda << " has support 2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << supp.k1
	   << ","
	   << supp.k2
	   << "]"
	   << endl;

      cout << "support intersection with first generator on level j0: ";
      bool inter = intersect_supports(basis, lambda, first_generator(&basis, basis.j0()), supp);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;
      
      if (lambda == last_wavelet(&basis, basis.j0())) break;
    }
#endif

#if 0
  cout << "- compute all intersecting wavelets:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0());; ++lambda)
    {
      cout << "  * for lambda=" << lambda << ":" << endl;
      typedef std::list<std::pair<Index, Basis::Support> > SupportList;
      SupportList nus;
      intersecting_wavelets(basis, lambda, basis.j0(), true, nus);
      for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	cout << "    nu=" << it->first 
	     << " with support intersection "
	     << "2^{-" << it->second.j << "}[" << it->second.k1 << "," << it->second.k2 << "]" << endl;
      }
      for (int level = basis.j0(); level <= basis.j0()+2; level++) {
	intersecting_wavelets(basis, lambda, level, false, nus);
	for (SupportList::const_iterator it(nus.begin()); it != nus.end(); ++it) {
	  cout << "    nu=" << it->first 
	       << " with support intersection "
	       << "2^{-" << it->second.j << "}[" << it->second.k1 << "," << it->second.k2 << "]" << endl;
	}
      }
      
      if (lambda == last_wavelet(&basis, basis.j0())) break;
    }
#endif  

#if 0
  cout << "- checking intersection of singular supports:" << endl;
  for (Index lambda = first_generator(&basis, basis.j0()+2);; ++lambda)
    {
      Support supp;
      support(basis, lambda, supp.k1, supp.k2);
      cout << "psi_lambda, lambda=" << lambda << " has the support 2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << supp.k1
	   << ","
	   << supp.k2
	   << "]"
	   << endl;
      
      Support supp_0;
      support(basis, first_generator(&basis, basis.j0()), supp_0.k1, supp_0.k2);
      cout << "* first generator on level j0 has the support 2^{-"
	   << basis.j0()
	   << "}["
	   << supp_0.k1
	   << ","
	   << supp_0.k2
	   << "]"
	   << endl;

      cout << "* support intersection with first generator on level j0         : ";
      bool inter = intersect_supports(basis, lambda, first_generator(&basis, basis.j0()), supp);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;

      cout << "* singular support intersection with first generator on level j0: ";
      inter = intersect_singular_support(basis, lambda, first_generator(&basis, basis.j0()), supp.j, supp.k1, supp.k2);
      if (inter)
	cout << "2^{-" << supp.j << "}[" << supp.k1 << "," << supp.k2 << "]" << endl;
      else
	cout << "none" << endl;

      if (lambda == last_wavelet(&basis, basis.j0()+2)) break;
    }
#endif

  return 0;
}

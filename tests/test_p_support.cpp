#include <iostream>

#include <interval/p_basis.h>
#include <interval/p_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of [P] generators and wavelets..." << endl;

  const int d = 2;
  const int dT = 2;

  typedef PBasis<d,dT> Basis;
  typedef Basis::Index Index;
  typedef Basis::Support Support;

//   Basis basis(0, 0); // no b.c.'s
//   Basis basis(1, 0); // complementary b.c. at x=0
//   Basis basis(0, 1); // complementary b.c. at x=1
  Basis basis(1, 1); // complementary b.c.'s

#if 1
  for (int level = basis.j0(); level <= basis.j0(); level++) {
    cout << "- computing the supports of all generators on level j=" << level << ":" << endl;
    
    Index lambda(first_generator(&basis, level));
    for (;; ++lambda) {
      int k1, k2;
      support(basis, lambda, k1, k2);
      cout << "  lambda=" << lambda << ", supp(psi_lambda)=2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << k1
	   << ","
	   << k2
	   << "]"
	   << endl;
      
      if (lambda == last_generator(&basis, level)) break;
    }
  }
#endif
  
#if 1
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

  return 0;
}

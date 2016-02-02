#include <iostream>

#include <Rd/cdf_basis.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Test RBasis support() routine..." << endl;

  const int d = 2;
  const int dt = 4;

  typedef CDFBasis<d,dt> Basis;
  typedef Basis::Index Index;
  typedef Basis::Support Support;
  Basis basis;

#if 1
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    cout << "- computing the supports of some generators on level j=" << level << ":" << endl;

    Index lambda(level,0,-10);
    for (;; ++lambda) {
      int k1, k2;
      basis.support(lambda, k1, k2);
      cout << "  lambda=" << lambda << ", supp(psi_lambda)=2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << k1
	   << ","
	   << k2
	   << "]"
	   << endl;
      
      if (lambda == Index(level,0,10)) break;
    }

    cout << "- computing the supports of some wavelets on level j=" << level << ":" << endl;

    lambda = Index(level,1,-15);
    for (;; ++lambda) {
      int k1, k2;
      basis.support(lambda, k1, k2);
      cout << "  lambda=" << lambda << ", supp(psi_lambda)=2^{-"
	   << lambda.j()+lambda.e()
	   << "}["
	   << k1
	   << ","
	   << k2
	   << "]"
	   << endl;
      
      if (lambda == Index(level,1,15)) break;
    }
  }
#endif

  return 0;
}

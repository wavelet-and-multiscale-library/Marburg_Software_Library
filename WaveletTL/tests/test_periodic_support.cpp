#include <iostream>

#include <Rd/cdf_basis.h>
#include <interval/periodic.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of periodic spline generators and wavelets..." << endl;

  const int d  = 2;
  const int dt = 2;

  typedef PeriodicBasis<CDFBasis<d,dt> > Basis;
  typedef Basis::Index Index;
  typedef Basis::Support Support;
  Basis basis;

#if 1
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    cout << "- computing the supports of all generators and wavelets on level j=" << level << ":" << endl;
    
    Index lambda(basis.first_generator(level));
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
      
      if (lambda == basis.last_wavelet(level)) break;
    }
  }
#endif
  
  return 0;
}

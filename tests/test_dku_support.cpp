#include <iostream>

#include <interval/dku_basis.h>
#include <interval/dku_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of DKU wavelets..." << endl;

  const int d = 2;
  const int dT = 2;

  typedef DKUBasis<d,dT> Basis;
  typedef Basis::Index Index;
  
  Basis basis;
  Index lambda(basis.firstGenerator(basis.j0()));
//   ++lambda;
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = basis.lastGenerator(basis.j0());
//   lambda = basis.firstWavelet(basis.j0()+1);
//   for (int i = 1; i <= 4; i++, ++lambda);
  lambda = basis.lastWavelet(basis.j0()+1);

  cout << "- point values at dyadic points of the " << lambda << " generator/wavelet:" << endl;
  basis.evaluate(lambda, true, 5).matlab_output(cout);

  int k1, k2;
  support(basis, lambda, k1, k2);
  cout << "- support of psi_lambda is 2^{-"
       << lambda.j()+lambda.e()
       << "}["
       << k1
       << ","
       << k2
       << "]"
       << endl;
  
  return 0;
}

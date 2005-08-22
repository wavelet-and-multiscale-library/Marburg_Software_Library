#include <iostream>

#include <interval/dku_basis.h>
#include <interval/dku_evaluate.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing point evaluation of DKU bases..." << endl;

  const int d = 2;
  const int dT = 4;

  typedef DKUBasis<d,dT> Basis;
  typedef Basis::Index Index;
  
//   Basis basis;
  Basis basis(true, true);
  Index lambda(basis.firstGenerator(basis.j0()));
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = basis.lastGenerator(basis.j0());
//   lambda = basis.firstWavelet(basis.j0());

  cout << "- point values at dyadic points of the " << lambda << " generator/wavelet:" << endl;
  basis.evaluate(lambda, true, 6).matlab_output(cout);

  cout << "- point values with evaluate() routine:" << endl;
  for (double x = 0.0; x <= 1.0; x += ldexp(1.0, -5))
    {
      cout << "x=" << x << ", psi(x)="
	   << evaluate(basis, 0, lambda, x) << endl;    
    }
  
  return 0;
}

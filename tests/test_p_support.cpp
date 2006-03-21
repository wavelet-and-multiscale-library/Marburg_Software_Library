#include <iostream>

#include <interval/p_basis.h>
#include <interval/p_support.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing support calculation of [P] generators and wavelets..." << endl;

  const int d = 2;
  const int dT = 4;

  typedef PBasis<d,dT> Basis;
  typedef Basis::Index Index;
  typedef Basis::Support Support;
  
//   Basis basis;
//   Index lambda(first_generator(&basis, basis.j0()));
//   ++lambda;
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = last_generator(&basis, basis.j0());
//   lambda = first_wavelet(&basis, basis.j0()+1);
//   for (int i = 1; i <= 4; i++, ++lambda);
//   lambda = last_wavelet(&basis, basis.j0()+1);

  return 0;
}

#include <iostream>

#include <algebra/infinite_vector.h>

#include <interval/ds_basis.h>
#include <Ldomain/ldomain_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the L-shaped domain wavelet index class ..." << endl;

  const int d  = 2;
  const int dT = 2;

  typedef DSBasis<d,dT> Basis1D;
  LDomainBasis<Basis1D> basis;

  LDomainIndex<Basis1D> index;
  cout << index << endl;

//   for (int i = 0; i < 64; i++) 
//     {
//       cout << index << endl;
//       ++index;
//     }
}

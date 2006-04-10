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
  typedef LDomainBasis<Basis1D> Basis;
  Basis basis;

  typedef Basis::Index Index;

  // the following output has to be checked manually...
  for (Index lambda(first_generator<Basis1D>(&basis, basis.j0()));; ++lambda) {
    cout << lambda << endl;
//     if (lambda.j() == 4) break;
//     if (lambda == last_generator<Basis1D>(&basis, basis.j0())) break;
//     if (lambda == first_wavelet<Basis1D>(&basis, basis.j0())) break;
    if (lambda == last_wavelet<Basis1D>(&basis, basis.j0()+1)) break;
  }
	 
}

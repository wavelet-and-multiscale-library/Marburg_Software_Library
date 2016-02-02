#include <iostream>

#include <algebra/infinite_vector.h>

#include <interval/ds_basis.h>
#include <interval/spline_basis.h>
#include <Ldomain/ldomain_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the L-shaped domain wavelet index class ..." << endl;

  const int d  = 2;
  const int dT = 2;

#if 1
  typedef DSBasis<d,dT,BernsteinSVD> Basis1D;
  Basis1D basis1d;
#else
  // this branch does not work at the moment!!!
  typedef SplineBasis<d,dT,DS_construction_bio5e> Basis1D;
  Basis1D basis1d;
#endif
  typedef LDomainBasis<Basis1D> Basis;
  Basis basis(basis1d);

  typedef Basis::Index Index;

#if 0
  // the following output has to be checked manually...
  for (Index lambda(first_generator<Basis1D>(&basis, basis.j0()));; ++lambda) {
    cout << lambda << endl;
//     if (lambda.j() == 4) break;
//     if (lambda == last_generator<Basis1D>(&basis, basis.j0())) break;
//     if (lambda == first_wavelet<Basis1D>(&basis, basis.j0())) break;
    if (lambda == last_wavelet<Basis1D>(&basis, basis.j0()+1)) break;
  }
#endif

#if 1
  cout << "- testing the number() routine:" << endl;
  int no = 0;
  for (Index lambda(basis.first_generator(basis.j0()));; ++lambda, ++no) {
    int number = lambda.number();
    cout << "lambda=" << lambda << " has the number " << number;
    if (number == no)
      cout << " (OK)" << endl;
    else
      cout << " (NO!!! I expect number=" << no << " here!)" << endl;
    if (lambda == basis.last_wavelet(basis.j0())) break;
  }
#endif
	 
}

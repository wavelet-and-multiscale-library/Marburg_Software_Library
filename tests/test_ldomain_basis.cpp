#include <iostream>

#include <algebra/infinite_vector.h>

#include <interval/ds_basis.h>
#include <Ldomain/ldomain_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the L-shaped domain..." << endl;

  const int d  = 2;
  const int dT = 2;

  typedef DSBasis<d,dT> Basis1D;
  typedef LDomainBasis<Basis1D> Basis;
  Basis basis;

  typedef Basis::Index Index;

#if 1
  cout << "- evaluating a primal generator..." << endl;
  Index lambda(first_generator<Basis1D>(&basis, basis.j0()));
  std::ofstream psistream("Ldomain_wavelet.m");
  matlab_output(psistream, evaluate<Basis1D>(basis, lambda, true, 6));
  psistream.close();
  cout << "  ...done, see file cube_wavelet.m!" << endl;
#endif

}

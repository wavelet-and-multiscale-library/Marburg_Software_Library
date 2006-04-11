#include <iostream>

#include <algebra/infinite_vector.h>

#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <Ldomain/ldomain_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the L-shaped domain..." << endl;

  const int d  = 2;
  const int dT = 2;

//   typedef DSBasis<d,dT> Basis1D; // remember to set partialSVD/BernsteinSVD biorthogonalization here!
  typedef PBasis<d,dT> Basis1D;
  typedef LDomainBasis<Basis1D> Basis;
  Basis basis;

  typedef Basis::Index Index;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- the default wavelet index w.r.t. the cube basis: " << Index(&basis) << endl;
  cout << "- first generator on the coarsest level: " << first_generator<Basis1D>(&basis, basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << last_generator<Basis1D>(&basis, basis.j0()) << endl;
  cout << "- first wavelet on the coarsest level: " << first_wavelet<Basis1D>(&basis, basis.j0()) << endl;
  cout << "- last wavelet on the coarsest level: " << last_wavelet<Basis1D>(&basis, basis.j0()) << endl;

#if 1
  Index lambda(first_generator<Basis1D>(&basis, basis.j0()));
  for (int i = 0; i < 273; i++, ++lambda);
  cout << "- evaluating a primal generator lambda=" << lambda << " ..." << endl;
  std::ofstream psistream("Ldomain_wavelet.m");
  matlab_output(psistream, evaluate<Basis1D>(basis, lambda, true, 6));
  psistream.close();
  cout << "  ...done, see file Ldomain_wavelet.m!" << endl;
#endif

}

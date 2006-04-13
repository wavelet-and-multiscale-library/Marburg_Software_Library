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

//   typedef DSBasis<d,dT> Basis1D; // remember to set partialSVD/BernsteinSVD biorthogonalization here! (why, comment this...)
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
//   for (; !(lambda.p() == 1); ++lambda);
//   for (; !(lambda.p() == 2); ++lambda);
//   for (; !(lambda.p() == 3); ++lambda);
  for (; !(lambda.p() == 4); ++lambda);
//   for (; lambda.e()[0] != 0 || lambda.e()[1] != 1; ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 1); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 1 && lambda.k()[1] == 7); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 2); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 2 && lambda.k()[1] == 7); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 4); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 4 && lambda.k()[1] == 6); ++lambda);
//   for (; !(lambda.e()[0] == 0 && lambda.e()[1] == 1 && lambda.p() == 4 && lambda.k()[1] == 7); ++lambda);

//   for (int i = 0; i < 155; i++, ++lambda); // one of the generators on patch 4
//   for (int i = 0; i < 330; i++, ++lambda); // one of the (0,1)-wavelets on patch 4
//   for (int i = 0; i < 334; i++, ++lambda);
  cout << "- evaluating a primal generator lambda=" << lambda << " ..." << endl;
  std::ofstream psistream("Ldomain_wavelet.m");
  matlab_output(psistream, evaluate<Basis1D>(basis, lambda, true, 6));
  psistream.close();
  cout << "  ...done, see file Ldomain_wavelet.m!" << endl;
#endif

}

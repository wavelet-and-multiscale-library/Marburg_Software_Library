#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algebra/infinite_vector.h>
#include <interval/periodic.h>
#include <ring/ring_basis.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;

  const double r0 = 0.5;
  const double r1 = 2.0;

  typedef RingBasis<d,dt,1,0> Basis;
  Basis basis(r0, r1);

  typedef Basis::Index Index;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- first generator on the coarsest level: " << basis.first_generator(basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << basis.last_generator(basis.j0()) << endl;
  cout << "- first wavelet on the coarsest level: " << basis.first_wavelet(basis.j0()) << endl;
  cout << "- last wavelet on the coarsest level: " << basis.last_wavelet(basis.j0()) << endl;

#if 0
  cout << "- testing iterator functionality:" << endl;
  for (Index lambda = basis.first_generator(basis.j0());; ++lambda) {
    cout << lambda << endl;
    if (lambda == basis.last_wavelet(basis.j0()+1)) break;
  }
#endif

#if 1
  const int resi = 6;
  Point<2> sphi, y;
  Matrix<double> gridx((1<<resi)+1), gridy((1<<resi)+1);
  for (int i = 0; i <= 1<<resi; i++) {
    sphi[0] = i/(double)(1<<resi);
    for (int j = 0; j <= 1<<resi; j++) {
      sphi[1] = j/(double)(1<<resi);
      basis.chart().map_point(sphi, y);
      gridx.set_entry(j, i, y[0]);
      gridy.set_entry(j, i, y[1]);
    }
  }
  Grid<2> grid(gridx, gridy);

  cout << "- evaluating a primal generator..." << endl;
  Index lambda(basis.first_generator(basis.j0()));
//   for (int i = 1; i <= 6; i++, ++lambda);
  cout << "  (lambda=" << lambda << " )" << endl;
  std::ofstream psistream("ring_wavelet.m");
  SampledMapping<2> psism = basis.evaluate(lambda, resi);
  if (false) {
    // plot in the parameter domain
    psism.matlab_output(psistream);
  } else {
    // plot in "world coordinates" on the ring
    SampledMapping<2> help(grid,psism.values());
    help.matlab_output(psistream);
  }
  psistream.close();
  cout << "  ...done, see file ring_wavelet.m!" << endl;
  
  cout << "- evaluating a primal wavelet:" << endl;
  lambda = ++basis.first_wavelet(basis.j0());
  cout << "  (lambda=" << lambda << " )" << endl;
  InfiniteVector<double,Index> coeffs;
  coeffs[lambda] = 1.0;
  std::ofstream psi2stream("ring_wavelet2.m");
  SampledMapping<2> psi2sm = basis.evaluate(coeffs, resi);
  if (false) {
    // plot in the parameter domain
    psi2sm.matlab_output(psi2stream);
  } else {
    // plot in "world coordinates" on the ring
    SampledMapping<2> help(grid,psi2sm.values());
    help.matlab_output(psi2stream);
  }
  psi2stream.close();
  cout << "  ...done, see file ring_wavelet2.m!" << endl;
#endif

#if 0
  cout << "- evaluating a whole lot of primal generators/wavelets..." << endl;
  Index mu(basis.first_generator(basis.j0()));
  for (int i = 1; i <= 64; i++, ++mu) {
    cout << "  * mu=" << mu << endl;
    ostringstream os;
    os << "ring_wavelet_" << mu << ".m";
    ofstream psistream(os.str().c_str());
    basis.evaluate(mu, 6).matlab_output(psistream);
    psistream.close();
  }
#endif
  
  return 0;
}

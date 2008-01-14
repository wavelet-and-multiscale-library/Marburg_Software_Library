#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <interval/periodic.h>
#include <ring/ring_basis.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;

  typedef RingBasis<d,dt,1,0> Basis;
  Basis basis;

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

#if 0
  cout << "- evaluating a primal generator..." << endl;
  Index lambda(basis.first_generator(basis.j0()));
//   for (int i = 1; i <= 6; i++, ++lambda);
  cout << "  (lambda=" << lambda << " )" << endl;
  std::ofstream psistream("ring_wavelet.m");
  basis.evaluate(lambda, 6).matlab_output(psistream);
  psistream.close();
  cout << "  ...done, see file ring_wavelet.m!" << endl;
  
  cout << "- evaluating a primal wavelet:" << endl;
  lambda = ++basis.first_wavelet(basis.j0());
  cout << "  (lambda=" << lambda << " )" << endl;
  std::ofstream psistream2("ring_wavelet2.m");
  basis.evaluate(lambda, 6).matlab_output(psistream2);
  psistream2.close();
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

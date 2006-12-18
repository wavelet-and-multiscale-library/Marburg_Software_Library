#include <iostream>
#include <fstream>

#include <Ldomain/ldomain_jl_basis.h>
#include <Ldomain/ldomain_jl_support.h>
#include <Ldomain/ldomain_jl_evaluate.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the JL basis on the L-shaped domain..." << endl;
  
  typedef LDomainJLBasis Basis;
  typedef Basis::Index Index;

  Basis basis;

  cout << "- j0=" << basis.j0() << endl;
  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- first generator on the coarsest level: " << first_generator(basis.j0()) << endl;
  cout << "- last generator on the coarsest level: " << last_generator(basis.j0()) << endl;
//   cout << "- first wavelet on the coarsest level: " << first_wavelet<Basis1D>(&basis, basis.j0()) << endl;
//   cout << "- last wavelet on the coarsest level: " << last_wavelet<Basis1D>(&basis, basis.j0()) << endl;

  int testlevel = 2;
  cout << "- first generator on the level " << testlevel << ": " << first_generator(testlevel) << endl;
  cout << "- last generator on the level " << testlevel << ": " << last_generator(testlevel) << endl;
//   cout << "- leftmost wavelet on the level " << testlevel << ": " << first_wavelet(testlevel) << endl;
//   cout << "- rightmost wavelet on the level " << testlevel << ": " << last_wavelet(testlevel) << endl;

#if 0
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    cout << "- iterate through all generators and wavelets on level j=" << level << ":" << endl;
    
//     Index lambda(first_generator(level));
//     Index lambda(first_wavelet(level));
    Index lambda(first_wavelet(level,Index::type_type(1,1)));
    for (;; ++lambda) {
      cout << lambda << endl;

      if (lambda == last_wavelet(level)) break;
//       if (lambda == last_generator(level)) break;
    }
  }
#endif

#if 1
//   Index lambda = first_generator(basis.j0());
  Index lambda = first_wavelet(basis.j0());
//   for (int i = 1; i <= 36; i++) ++lambda;

  cout << "- evaluating psi_{" << lambda << "}..." << endl;
  std::ofstream psistream("Ldomain_wavelet.m");
  matlab_output(psistream, evaluate(basis, lambda, 6));
  psistream.close();
  cout << "  ...done, see file Ldomain_wavelet.m!" << endl;
#endif

  return 0;
}
 

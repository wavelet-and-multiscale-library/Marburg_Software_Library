#include <iostream>

#include <Ldomain/ldomain_jl_basis.h>

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



  return 0;
}

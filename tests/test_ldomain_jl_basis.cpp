#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>

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
  for (int level = basis.j0(); level <= basis.j0(); level++) {
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
  const int j1 = basis.j0()+1;
  // prepare set of valid wavelet indices
  set<Index> valid;
  for (Index lambda = first_generator(j1); lambda <= last_generator(j1); ++lambda) {
    valid.insert(lambda);
  }
  cout << "- a set of 'valid' wavelet indices:" << endl;
  for (set<Index>::const_iterator it = valid.begin(); it != valid.end(); ++it)
    cout << *it << endl;

  set<Index> tested;

  for (int c0 = 0; c0 <= 1; c0++)
    for (int c1 = 0; c1 <= 1; c1++)
      for (int k0 = -(1<<j1); k0 <= (1<<j1); k0++)
	for (int k1 = -(1<<j1); k1 <= (1<<j1); k1++) {
	  if (index_is_valid(j1,0,0,c0,c1,k0,k1))
	    tested.insert(Index(j1,Index::type_type(0,0),Index::component_type(c0,c1),Index::translation_type(k0,k1)));
	}
  cout << "- a set of 'tested' wavelet indices:" << endl;
  for (set<Index>::const_iterator it = tested.begin(); it != tested.end(); ++it)
    cout << *it << endl;

  set<Index> tested_but_invalid;
  set_difference(tested.begin(), tested.end(),
		 valid.begin(), valid.end(),
		 inserter(tested_but_invalid, tested_but_invalid.begin()));
  cout << "- 'tested' but really invalid wavelet indices:" << endl;
  for (set<Index>::const_iterator it = tested_but_invalid.begin(); it != tested_but_invalid.end(); ++it)
    cout << *it << endl;

  set<Index> valid_but_not_tested;
  set_difference(valid.begin(), valid.end(),
		 tested.begin(), tested.end(),
		 inserter(valid_but_not_tested, valid_but_not_tested.begin()));
  cout << "- valid but wrongly tested wavelet indices:" << endl;
  for (set<Index>::const_iterator it = valid_but_not_tested.begin(); it != valid_but_not_tested.end(); ++it)
    cout << *it << endl;

#endif  

#if 0
  Index lambda = first_generator(basis.j0());
//   Index lambda = first_wavelet(basis.j0());
//   for (int i = 1; i <= 36; i++) ++lambda;
  for (int i = 1; i <= 24; i++) ++lambda;

  cout << "- evaluating psi_{" << lambda << "}..." << endl;
  std::ofstream psistream("Ldomain_wavelet.m");
  matlab_output(psistream, evaluate(basis, lambda, 6));
  psistream.close();
  cout << "  ...done, see file Ldomain_wavelet.m!" << endl;
#endif

  return 0;
}
 

#include <iostream>
#include <interval/jl_basis.h>
#include <interval/jl_evaluate.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases from [JL] ..." << endl;

  typedef JLBasis Basis;
  typedef Basis::Index Index;

  Basis basis;

  cout << "- the default wavelet index: " << Index() << endl;
  cout << "- leftmost generator on the coarsest level: " << first_generator(basis.j0()) << endl;
  cout << "- rightmost generator on the coarsest level: " << last_generator(basis.j0()) << endl;
  cout << "- leftmost wavelet on the coarsest level: " << first_wavelet(basis.j0()) << endl;
  cout << "- rightmost wavelet on the coarsest level: " << last_wavelet(basis.j0()) << endl;

  int testlevel = 2;
  cout << "- leftmost generator on the level " << testlevel << ": " << first_generator(testlevel) << endl;
  cout << "- rightmost generator on the level " << testlevel << ": " << last_generator(testlevel) << endl;
  cout << "- leftmost wavelet on the level " << testlevel << ": " << first_wavelet(testlevel) << endl;
  cout << "- rightmost wavelet on the level " << testlevel << ": " << last_wavelet(testlevel) << endl;


#if 0
  for (int level = basis.j0(); level <= basis.j0()+1; level++) {
    cout << "- iterate through all generators and wavelets on level j=" << level << ":" << endl;
    
    Index lambda(first_generator(level));
    for (;; ++lambda) {
      cout << lambda << endl;

      if (lambda == last_wavelet(level)) break;
//       if (lambda == last_generator(level)) break;
    }
  }
#endif
  
#if 0
//   Index mu(1,0,1,0);
//   Index mu(1,1,1,0);
//   Index mu(first_wavelet(basis.j0()));
  Index mu(first_generator(basis.j0())); ++mu;
  cout << "* a wavelet coefficient mu=" << mu << endl;
  cout << "  - reconstruct_1() yields" << endl;
  InfiniteVector<double,Index> gcoeffs;
  basis.reconstruct_1(mu,mu.j()+1,gcoeffs);
  cout << gcoeffs;
#endif

#if 1
  cout << "- evaluating some primal generators:" << endl;
  for (Index lambda(first_generator(basis.j0()));; ++lambda) {
    cout << lambda << endl;
    evaluate(basis, lambda, 6).matlab_output(cout);
    if (lambda == last_generator(basis.j0())) break;
  }
  
  cout << "- evaluating some primal wavelets:" << endl;
  for (Index lambda = first_wavelet(basis.j0());; ++lambda) {
    cout << lambda << endl;
    evaluate(basis, lambda, 6).matlab_output(cout);
    if (lambda == last_wavelet(basis.j0()+1)) break;
  }
#endif

  return 0;
}

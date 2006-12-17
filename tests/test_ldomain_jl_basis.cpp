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


  return 0;
}

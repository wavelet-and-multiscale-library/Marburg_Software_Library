#include <iostream>
#include <fstream>
#include <interval/spline_basis_data.h>
#include <galerkin/full_laplacian.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing FullLaplacian ..." << endl;

  SplineBasisData<2,2> sd22("P","",1,1,0,0); // PBasis, complementary b.c.'s
  FullLaplacian<2,2> delta(sd22);

  return 0;
}

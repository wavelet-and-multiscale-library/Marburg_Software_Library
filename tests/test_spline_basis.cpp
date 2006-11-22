#include <iostream>
#include <fstream>
#include <interval/spline_basis_data.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing setup of SplineBasisData objects ..." << endl;

  const unsigned int d = 3;
  const unsigned int dT = 3;

  SplineBasisData<d,dT> sd("P","",1,1,0,0); // PBasis, d=2, dT=2, complementary b.c.'s
  sd.check();

  return 0;
}

#include <iostream>
#include <interval/dku.h>
#include <utils/array1d.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DKU bases..." << endl;

  const int d = 3;
  const int dT = 5;
  DKUBasis<d, dT> basis;
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;

  return 0;
}

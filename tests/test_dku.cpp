#include <iostream>
#include <interval/dku.h>
#include <utils/array1d.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DKU bases..." << endl;

  const int d = 2;
  const int dT = 2;
  DKUBasis<d, dT> basis(none);
  cout << "- the (" << d << "," << dT << ") basis has j0=" << basis.j0() << endl;

  return 0;
}

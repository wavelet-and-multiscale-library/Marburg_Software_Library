#include <iostream>
#include <interval/dku.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DKU bases..." << endl;

  DKUBasis<2, 2> basis;
  cout << "- the (2,2) basis has j0=" << basis.j0() << endl;

  return 0;
}

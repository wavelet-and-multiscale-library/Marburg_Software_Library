#include <iostream>
#include <interval/dku.h>
#include <utils/array1d.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the DKU bases..." << endl;

  Array1D<int> bc(2);
  bc[0] = bc[1] = 0; // free boundary conditions
  DKUBasis<2, 2> basis(bc, bc);
  cout << "- the (2,2) basis has j0=" << basis.j0() << endl;

  return 0;
}

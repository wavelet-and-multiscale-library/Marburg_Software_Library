#include <iostream>
#include <Rd/haar_mask.h>
#include <Rd/dm_mask.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing the Dahmen/Micchelli mask ..." << endl;

  DMMask<HaarMask, HaarMask> dm;
  cout << dm << endl;

  return 0;
}

#include <iostream>
#include <Rd/haar_mask.h>
#include <interval/periodic.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing periodic wavelet bases ..." << endl;

  PeriodicBasis<HaarMask> mask;

  return 0;
}

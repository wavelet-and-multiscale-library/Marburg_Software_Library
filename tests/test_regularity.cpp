#include <iostream>

#include <utils/multiindex.h>

#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/daubechies_mask.h>
#include <Rd/refinable.h>
#include <Rd/regularity.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;

int main()
{
  cout << "Computing the (critical) Sobolev regularity of some univariate refinable functions..." << endl;

  cout << "- Sobolev regularity of the Haar function: " << endl;
  cout << Sobolev_regularity<HaarMask>() << endl;

  return 0;
}

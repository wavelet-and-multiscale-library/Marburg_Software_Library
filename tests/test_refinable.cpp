#include <iostream>
#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/refinable.h>

using namespace std;
using namespace WaveletTL;

int main()
{
  cout << "Testing refinable functions..." << endl;

  cout << "- evaluate a Haar function:" << endl;
  RefinableFunction<HaarMask> phi;
  cout << "  + primal mask:" << endl;
  for (RefinableFunction<HaarMask>::const_iterator it(phi.begin()); it != phi.end(); it++)
    cout << "k=" << it.index() << ": " << *it << endl;
  phi.evaluate<0u>(0, 0, -2, 2, 1).matlab_output(cout);

  cout << "- evaluate a primal CDF<2> function:" << endl;
  RefinableFunction<CDFMask_primal<2> > phi2;
  for (RefinableFunction<CDFMask_primal<2> >::const_iterator it(phi2.begin()); it != phi2.end(); it++)
    cout << "k=" << it.index() << ": " << *it << endl;
  phi2.evaluate<0u>(0, 0, -1, 2, 2).matlab_output(cout);

  cout << "- evaluate first derivative of this function:" << endl;
  phi2.evaluate<1u>(0, 0, -1, 2, 2).matlab_output(cout);

  cout << "- CDF<3>, function values and first derivative:" << endl;
  RefinableFunction<CDFMask_primal<3> > phi3;
  phi3.evaluate<0u>(0, 0, -1, 2, 2).matlab_output(cout);
  phi3.evaluate<1u>(0, 0, -1, 2, 2).matlab_output(cout);

  cout << "- hat function and derivative on level 1:" << endl;
  phi2.evaluate<0u>(1, 0, -1, 2, 2).matlab_output(cout);
  phi2.evaluate<1u>(1, 0, -1, 2, 2).matlab_output(cout);

  cout << "- calculating moments of a CDF<2> function:" << endl;
  for (unsigned int n(0); n <= 9; n++)
    {
      cout << "    n=" << n << ": M_n=" << phi2.moment(n)
  	   << " (exact value: ";
      if (n % 2 == 0)
  	cout << 2.0/((n+1)*(n+2)); // exact moments, cf. [BDK]
      else
  	cout << "0";
      cout << ")" << endl;
    }

  cout << "- evaluate the CDF<2> function with the [DM] version of evaluate():" << endl
       << phi2.evaluate(0) << endl;

  cout << "- evaluate the CDF<3> function with the [DM] version of evaluate():" << endl
       << phi3.evaluate(0) << endl;

  return 0;
}

#include <iostream>

#include <utils/multiindex.h>

#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/refinable.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;

int main()
{
  cout << "Testing refinable functions..." << endl;

  cout << "- evaluate a Haar function:" << endl;
  RefinableFunction<HaarMask> phi;
  cout << "  + primal mask:" << endl;
  for (RefinableFunction<HaarMask>::const_iterator it(phi.begin()); it != phi.end(); it++)
    cout << "k=" << it.index() << ": " << *it << endl;
  cout << "  + integer values:" << endl;
  cout << phi.evaluate();
  cout << "  + values at a finer resolution:" << endl;
  cout << phi.evaluate(2);

  cout << "- evaluate a primal CDF<2> function:" << endl;
  RefinableFunction<CDFMask_primal<2> > phi2;
  for (RefinableFunction<CDFMask_primal<2> >::const_iterator it(phi2.begin()); it != phi2.end(); it++)
    cout << "k=" << it.index() << ": " << *it << endl;
  cout << "  + integer values:" << endl;
  cout << phi2.evaluate();
  cout << "  + values at a finer resolution:" << endl;
  cout << phi2.evaluate(2);

  cout << "- CDF<3>, function values and first derivative:" << endl;
  RefinableFunction<CDFMask_primal<3> > phi3;
  cout << phi3.evaluate(2);
  cout << phi3.evaluate(MultiIndex<unsigned int, 1>(1), 2);

//   cout << "- calculating moments of a CDF<2> function:" << endl;
//   for (unsigned int n(0); n <= 9; n++)
//     {
//       cout << "    n=" << n << ": M_n=" << phi2.moment(n)
//   	   << " (exact value: ";
//       if (n % 2 == 0)
//   	cout << 2.0/((n+1)*(n+2)); // exact moments, cf. [BDK]
//       else
//   	cout << "0";
//       cout << ")" << endl;
//     }

//   cout << "- evaluate the CDF<2> function with the [DM] version of evaluate():" << endl
//        << phi2.evaluate(0) << endl;

//   cout << "- evaluate the CDF<3> function with the [DM] version of evaluate():" << endl
//        << phi3.evaluate(0) << endl;

  return 0;
}

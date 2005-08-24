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

  cout << "- dual mask of a CDF<2,2> function:" << endl;
  RefinableFunction<CDFMask_dual<2,2> > phi2T;
  for (RefinableFunction<CDFMask_dual<2,2> >::const_iterator it(phi2T.begin()); it != phi2T.end(); it++)
    cout << "k=" << it.index() << ": " << *it << endl;  

  cout << "- CDF<3>, function values and first derivative, clipped to [0,1]:" << endl;
  RefinableFunction<CDFMask_primal<3> > phi3;
  cout << phi3.evaluate(2);
  cout << phi3.evaluate(MultiIndex<unsigned int, 1>(1), 2);

  cout << "- the Haar function on a higher level:" << endl;
  phi.evaluate(MultiIndex<unsigned int, 1>(),
	       1, MultiIndex<int, 1>(1),
	       MultiIndex<int, 1>(0), MultiIndex<int, 1>(1),
	       2).matlab_output(cout);

  cout << "- the CDF<3> function phi_{0,0}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<unsigned int, 1>(),
		0, MultiIndex<int, 1>(0),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);
  cout << "- the CDF<3> function phi_{1,-1}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<unsigned int, 1>(),
		1, MultiIndex<int, 1>(-1),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);
  cout << "- first derivative of the CDF<3> function phi_{0,0}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<unsigned int, 1>(1),
		0, MultiIndex<int, 1>(0),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);
  cout << "- first derivative of the CDF<3> function phi_{1,1}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<unsigned int, 1>(1),
		1, MultiIndex<int, 1>(1),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);

  cout << "- calculating moments of a CDF<2> function:" << endl;
  for (unsigned int n(0); n <= 9; n++)
    {
      cout << "    n=" << n << ": M_n="
	   << phi2.moment(MultiIndex<unsigned int, 1>(n))
   	   << " (exact value: ";
      if (n % 2 == 0)
   	cout << 2.0/((n+1)*(n+2)); // exact moments, cf. [BDK]
      else
   	cout << "0";
      cout << ")" << endl;
    }
  
  return 0;
}

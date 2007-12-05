#include <iostream>

#include <utils/multiindex.h>

#include <Rd/haar_mask.h>
#include <Rd/cdf_mask.h>
#include <Rd/daubechies_mask.h>
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
  for (int k = phi.abegin(); k <= phi.aend(); k++)
    cout << "k=" << k << ": " << phi.a(k) << endl;
  cout << "  + integer values:" << endl;
  cout << phi.evaluate();
  cout << "  + values at a finer resolution:" << endl;
  cout << phi.evaluate(2);
  cout << "  + evaluate() with 0-th derivative at a finer resolution:" << endl;
  cout << phi.evaluate(0, 2);
  cout << "  + evaluate() with 0-th derivative, j=0, k=0 at a finer resolution, clipped to [0,1]:" << endl;
  phi.evaluate(0,
 	       0, 0,
 	       0, 1,
 	       2).matlab_output(cout);
  cout << "  + evaluate() with 0-th derivative, j=1, k=1 at a finer resolution, clipped to [0,1]:" << endl;
  phi.evaluate(0,
 	       1, 1,
 	       0, 1,
 	       2).matlab_output(cout);
  
#if 1
  cout << "- evaluate a primal CDF<2> function:" << endl;
  RefinableFunction<CDFRefinementMask_primal<2> > phi2;
  for (int k = phi2.abegin(); k <= phi2.aend(); k++)
    cout << "k=" << k << ": " << phi2.a(k) << endl;
  cout << "  + integer values:" << endl;
  cout << phi2.evaluate();
  cout << "  + values at a finer resolution:" << endl;
  cout << phi2.evaluate(2);  
  cout << "  + evaluate() with 0-th derivative at a finer resolution:" << endl;
  cout << phi2.evaluate(0, 2);
#endif

#if 1
  cout << "- dual mask of a CDF<2,2> function:" << endl;
  RefinableFunction<CDFRefinementMask_dual<2,2> > phi2T;
  for (int k = phi2T.abegin(); k <= phi2T.aend(); k++)
    cout << "k=" << k << ": " << phi2T.a(k) << endl;  

  cout << "- dual mask of a CDF<2,4> function:" << endl;
  RefinableFunction<CDFRefinementMask_dual<2,4> > phi24T;
  for (int k = phi24T.abegin(); k <= phi24T.aend(); k++)
    cout << "k=" << k << ": " << phi24T.a(k) << endl;  
#endif

#if 1
  cout << "- primal CDF<3> mask:" << endl;
  RefinableFunction<CDFRefinementMask_primal<3> > phi3;
  for (int k = phi3.abegin(); k <= phi3.aend(); k++)
    cout << "k=" << k << ": " << phi3.a(k) << endl;  
  cout << "  + integer values:" << endl;
  cout << phi3.evaluate();
  cout << "  + values at a finer resolution:" << endl;
  cout << phi3.evaluate(2);  
  cout << "  + evaluate() with 0-th derivative at a finer resolution:" << endl;
  cout << phi3.evaluate(0, 2);
  cout << "  + evaluate() with 1st derivative at a finer resolution:" << endl;
  cout << phi3.evaluate(1, 2);
#endif

#if 1
  cout << "- dual CDF<3,3> mask:" << endl;
  RefinableFunction<CDFRefinementMask_dual<3,3> > phi3T0;
  for (int k = phi3T0.abegin(); k <= phi3T0.aend(); k++)
    cout << "k=" << k << ": " << phi3T0.a(k) << endl;  

  cout << "- dual CDF<3,5> mask:" << endl;
  RefinableFunction<CDFRefinementMask_dual<3,5> > phi3T;
  for (int k = phi3T.abegin(); k <= phi3T.aend(); k++)
    cout << "k=" << k << ": " << phi3T.a(k) << endl;  
#endif

#if 0
  cout << "- CDF<3>, function values and first derivative, clipped to [0,1]:" << endl;
  cout << phi3.evaluate(2);
  cout << phi3.evaluate(MultiIndex<int, 1>(1), 2);

  cout << "- the Haar function on a higher level:" << endl;
  phi.evaluate(MultiIndex<int, 1>(),
	       1, MultiIndex<int, 1>(1),
	       MultiIndex<int, 1>(0), MultiIndex<int, 1>(1),
	       2).matlab_output(cout);
#endif

#if 0
  cout << "- the CDF<3> function phi_{0,0}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<int, 1>(),
		0, MultiIndex<int, 1>(0),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);
  cout << "- the CDF<3> function phi_{1,-1}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<int, 1>(),
		1, MultiIndex<int, 1>(-1),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);
  cout << "- first derivative of the CDF<3> function phi_{0,0}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<int, 1>(1),
		0, MultiIndex<int, 1>(0),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);
  cout << "- first derivative of the CDF<3> function phi_{1,1}, clipped to [-2,2]:" << endl;
  phi3.evaluate(MultiIndex<int, 1>(1),
		1, MultiIndex<int, 1>(1),
		MultiIndex<int, 1>(-2), MultiIndex<int, 1>(2),
		2).matlab_output(cout);
#endif

#if 0
  cout << "- calculating moments of a CDF<2> function:" << endl;
  for (unsigned int n(0); n <= 9; n++)
    {
      cout << "    n=" << n << ": M_n="
	   << phi2.moment(MultiIndex<int, 1>(n))
   	   << " (exact value: ";
      if (n % 2 == 0)
   	cout << 2.0/((n+1)*(n+2)); // exact moments, cf. [BDK]
      else
   	cout << "0";
      cout << ")" << endl;
    }
#endif
  
#if 0
  cout << "- Daubechies<2> mask:" << endl;
  RefinableFunction<DaubechiesMask<2> > dau2;
  for (RefinableFunction<DaubechiesMask<2> >::const_iterator it(dau2.begin()); it != dau2.end(); it++)
    cout << "k=" << it.index() << ": " << *it << endl;  
#endif

  return 0;
}

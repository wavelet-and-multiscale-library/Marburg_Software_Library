#include <iostream>

#include <geometry/sampled_mapping.h>
#include <Rd/cdf_basis.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Test CDFBasis class..." << endl;

  const int d = 3;
  const int dt = 3;

  CDFBasis<d,dt> basis;

  cout << "A CDF basis with d=" << d << " and dt=" << dt
       << " has the following data:" << endl;
  cout << "  + primal mask:" << endl;
  cout << basis.a() << endl;
  cout << "  + dual mask:" << endl;
  cout << basis.aT() << endl;

  cout << "  + wavelet coefficients:" << endl;
  for (int k = 1-CDFBasis<d,dt>::dual_mask::end();
       k <= 1-CDFBasis<d,dt>::dual_mask::begin();
       k++) {
    cout << "b(" << k << ")=" << basis.b(k) << endl;
  }

#if 1
  cout << "point values:" << endl;

  Grid<1> grid(-4.0, 4.0, 512);
  Array1D<double> values;
  basis.evaluate(0, RIndex(3,1,0), grid.points(), values);
  SampledMapping<1>(grid,values).matlab_output(cout);
#endif

  return 0;
}

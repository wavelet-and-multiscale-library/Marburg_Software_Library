#include <iostream>

#include <geometry/sampled_mapping.h>
#include <Rd/cdf_basis.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Test CDFBasis class..." << endl;

  const int d = 2;
  const int dt = 4;

  CDFBasis<d,dt> basis;

  Grid<1> grid(-4.0, 4.0, 512);
  Array1D<double> values(grid.size());
  for (unsigned int k = 0; k < values.size(); k++)
    values[k] = basis.evaluate(0, RIndex(3,1,0), grid.points()[k]);
  SampledMapping<1>(grid,values).matlab_output(cout);

  return 0;
}

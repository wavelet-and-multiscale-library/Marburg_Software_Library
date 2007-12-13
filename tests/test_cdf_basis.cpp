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
  const int dt = 2;

  CDFBasis<d,dt> basis;

  Grid<1> grid(-4.0, 4.0, 500);
  Array1D<double> values(grid.size());
  for (unsigned int k = 0; k < values.size(); k++)
    values[k] = basis.evaluate(1, RIndex(2,1,-1), grid.points()[k]);
  SampledMapping<1>(grid,values).matlab_output(cout);

  return 0;
}

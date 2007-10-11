#include <iostream>
#include <numerics/up_function.h>
#include <geometry/sampled_mapping.h>

using namespace std;
using namespace MathTL;

int main()
{
  for (unsigned int j = 0; j <= 10; j++) {
    ApproximateUpFunction uj(j);

    cout << "point values of u_" << j << " on a dyadic subgrid of [-1,1]:" << endl;
    SampledMapping<1> g(Grid<1>(-1.0, 1.0, 10), uj);
    g.matlab_output(cout);
  }

  return 0;
}

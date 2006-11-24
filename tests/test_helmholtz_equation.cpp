#include <iostream>
#include <interval/spline_basis.h>
#include <galerkin/helmholtz_equation.h>
#include <utils/function.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing HelmholtzEquation ..." << endl;

  const unsigned int d = 2;
  const unsigned int dT = 2;

  HelmholtzEquation1D<d,dT> helmholtz(0);

  return 0;
}

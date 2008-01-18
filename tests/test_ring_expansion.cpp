#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algebra/infinite_vector.h>
#include <algebra/matrix.h>
#include <geometry/grid.h>
#include <interval/periodic.h>
#include <utils/function.h>
#include <ring/ring_basis.h>
#include "ring_functions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;

  const double r0 = 0.5;
  const double r1 = 2.0;

  typedef RingBasis<d,dt,1,1> Basis;
  Basis basis(r0, r1);

  typedef Basis::Index Index;

  const int j0 = basis.j0();
  const int jmax = j0;

  Function<2>* u = 0;
  const int testcase = 3;
  switch(testcase) {
  case 1:
    u = new RingFunction1(r0, r1);
    break;
  case 2:
    u = new RingFunction2(r0, r1);
    break;
  case 3:
    u = new RingFunction3(r0, r1);
    break;
  case 4:
    u = new RingFunction4();
    break;
  default:
    u = new RingFunction1(r0, r1);
  }

  const int resi = 6;
  Point<2> sphi, y;
  Matrix<double> gridx((1<<resi)+1), gridy((1<<resi)+1);
  for (int i = 0; i <= 1<<resi; i++) {
    sphi[0] = i/(double)(1<<resi);
    for (int j = 0; j <= 1<<resi; j++) {
      sphi[1] = j/(double)(1<<resi);
      basis.chart().map_point(sphi, y);
      gridx.set_entry(j, i, y[0]);
      gridy.set_entry(j, i, y[1]);
    }
  }
  Grid<2> grid(gridx, gridy);

#if 1
  {
    // plot the test function
    cout << "- plotting the test function..." << endl;
    SampledMapping<2> rf(grid, *u);
    std::ofstream rfstream("test_function.m");
    rf.matlab_output(rfstream);
    rfstream.close();
    cout << "  ...done, see file test_function.m!" << endl;
  }
#endif
  
  InfiniteVector<double,Index> coeffs;
  
#if 1
  basis.expand(u, true, jmax, coeffs);
  cout << "- integrals of the test function u against all wavelets up to level jmax="
       << jmax << ":" << endl
       << coeffs << endl;
#endif

#if 1
  basis.expand(u, false, jmax, coeffs);
  cout << "- approximate expansion coefficients of the test function in the primal basis:" << endl
       << coeffs << endl;
  if (true) {
    cout << "- compute L_infty error on a subgrid..." << endl;
    SampledMapping<2> rf(grid, *u);
    SampledMapping<2> rf0(basis.evaluate(coeffs, resi));
    cout << "  ... done, result: "
	 << maximum_norm(rf.values()-rf0.values())
	 << endl;

    cout << "- plotting the approximate expansion..." << endl;
    SampledMapping<2> rf1(grid,rf0.values());
    std::ofstream rf1stream("expansion.m");
    rf1.matlab_output(rf1stream);
    rf1stream.close();
    cout << "  ...done, see file expansion.m!" << endl;

    cout << "- plotting the pointwise error..." << endl;
    SampledMapping<2> rferr(grid,rf.values()-rf0.values());
    std::ofstream rferrstream("exp_error.m");
    rferr.matlab_output(rferrstream);
    rferrstream.close();
    cout << "  ...done, see file exp_error.m!" << endl;
    
  }
#endif
  
  if (u) delete u;
  
  return 0;
}

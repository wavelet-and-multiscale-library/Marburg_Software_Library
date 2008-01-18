#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

#include <algebra/sparse_matrix.h>
#include <numerics/iteratsolv.h>
#include <ring/ring_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/ring_laplacian.h>
#include "ring_functions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  cout << "Testing discretizations of the Laplacian with wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;
  const int s0 = 1;
  const int s1 = 1;

  const double r0 = 0.5;
  const double r1 = 2.0;

  typedef RingBasis<d,dt,s0,s1> Basis;
  Basis basis(r0, r1);

  typedef Basis::Index Index;

  const int j0 = basis.j0();
  const int jmax = j0+1;

  Function<2> *u = 0, *f = 0;
  const int testcase = 4;
  switch(testcase) {
  case 4:
    u = new RingFunction4(r0, r1);
    f = new RingFunction4_RHS(r0, r1);
    break;
  default:
    u = new RingFunction4(r0, r1);
    f = new RingFunction4_RHS(r0, r1);
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

#if 0
  {
    // plot solution and right-hand side
    cout << "- plotting the exact solution..." << endl;
    SampledMapping<2> u_sm(grid, *u);
    std::ofstream u_stream("solution.m");
    u_sm.matlab_output(u_stream);
    u_stream.close();
    cout << "  ...done, see file solution.m!" << endl;
    cout << "- plotting the right-hand side..." << endl;
    SampledMapping<2> f_sm(grid, *f);
    std::ofstream f_stream("rhs.m");
    f_sm.matlab_output(f_stream);
    f_stream.close();
    cout << "  ...done, see file rhs.m!" << endl;
  }
#endif

  typedef RingLaplacian<d,dt,s0,s1> Problem;
  Problem problem(basis, InfiniteVector<double,Index>());



  if (u) delete u;
  if (f) delete f;

  return 0;
}

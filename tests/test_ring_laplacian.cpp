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

  const int d  = 3;
  const int dt = 3;
  const int s0 = 1;
  const int s1 = 1;

  const double r0 = 0.5;
  const double r1 = 2.0;

  typedef RingBasis<d,dt,s0,s1> Basis;
  Basis basis(r0, r1);

  typedef Basis::Index Index;

  const int j0 = basis.j0();
  const int jmax = j0;

  Function<2> *uexact = 0, *f = 0;
  const int testcase = 4;
  switch(testcase) {
  case 4:
    uexact = new RingFunction4(r0, r1);
    f = new RingFunction4_RHS(r0, r1);
    break;
  default:
    uexact = new RingFunction4(r0, r1);
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

#if 1
  {
    // plot solution and right-hand side
    cout << "- plotting the exact solution..." << endl;
    SampledMapping<2> u_sm(grid, *uexact);
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

  cout << "- expand right-hand side..." << endl;
  InfiniteVector<double,Index> fcoeffs;
  basis.expand(f, true, jmax, fcoeffs);
  problem.set_rhs(fcoeffs);
  cout << "  ... done!" << endl;

  cout << "- set up index set of active wavelets..." << endl;
  set<Index> Lambda;
  for (Index lambda = basis.first_generator(basis.j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == basis.last_wavelet(jmax)) break;
  }
  cout << "  ... done!" << endl;

  cout << "- set up stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;

  tstart = clock();
  SparseMatrix<double> A;
  setup_stiffness_matrix(problem, Lambda, A);
  tend = clock();

  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- stiffness matrix A=" << endl << A << endl;

  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(problem, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- right hand side: " << b << endl;
  
  Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
  unsigned int iterations;
  CG(A, b, x, 1e-15, 200, iterations);
  
  //   cout << "- solution coefficients: " << x;
  cout << "- residual (infinity) norm: ";
  A.apply(x, err);
  err -= b;
  cout << linfty_norm(err) << endl;
   
#if 1
  {
    InfiniteVector<double,Index> ucoeffs;
    unsigned int i = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
      ucoeffs.set_coefficient(*it, x[i]);
    ucoeffs.scale(&problem, -1);

    cout << "- compute L_infty error on a subgrid..." << endl;
    SampledMapping<2> u_sm(grid, *uexact);
    SampledMapping<2> u0_sm(basis.evaluate(ucoeffs, resi));
    cout << "  ... done, result: "
	 << maximum_norm(u_sm.values()-u0_sm.values())
	 << endl;
    
    cout << "- plot point values of the solution:" << endl;
    SampledMapping<2> u0_sm2(grid,u0_sm.values());
    std::ofstream u0_stream("u0.m");
    u0_sm2.matlab_output(u0_stream);
    u0_stream.close();
    cout << "  ...done, see file u0.m!" << endl;

    cout << "- plot the pointwise error..." << endl;
    SampledMapping<2> u0_err(grid,u_sm.values()-u0_sm.values());
    std::ofstream u0_err_stream("u0_error.m");
    u0_err.matlab_output(u0_err_stream);
    u0_err_stream.close();
    cout << "  ...done, see file u0_error.m!" << endl;
  }
#endif


  if (uexact) delete uexact;
  if (f) delete f;

  return 0;
}

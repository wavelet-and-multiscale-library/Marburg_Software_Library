#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

#include <algebra/sparse_matrix.h>
#include <numerics/iteratsolv.h>
#include <ring/ring_basis.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/ring_gramian.h>
#include "ring_functions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  cout << "Testing Gramian matrices of wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;
  const int s0 = 1;
  const int s1 = 1;
//   const int J0 = SplineBasisData_j0<d,dt,P_construction,s0,s1,0,0>::j0;

  const double r0 = 0.5;
  const double r1 = 2.0;

  typedef RingBasis<d,dt,s0,s1> Basis;
  Basis basis(r0, r1);

  typedef Basis::Index Index;

  const int j0 = basis.j0();
  const int jmax = j0+3;

  cout << "j0=" << j0 << ", jmax=" << jmax << endl;

  Function<2>* u = 0;
  const int testcase = 4;
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
    u = new RingFunction4(r0, r1);
    break;
  case 5:
    u = new RingFunction4_RHS(r0, r1);
    break;
  default:
    u = new RingFunction1(r0, r1);
  }

  typedef RingGramian<d,dt,s0,s1> Problem;
  Problem problem(basis, InfiniteVector<double,Index>());

  cout << "- expand right-hand side..." << endl;
  InfiniteVector<double,Index> fcoeffs;
  basis.expand(u, true, jmax, fcoeffs);
  problem.set_rhs(fcoeffs);
  cout << "  ... done!" << endl;

  cout << "- set up index set of active wavelets..." << endl;
  set<Index> Lambda;
  for (Index lambda = basis.first_generator(basis.j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == basis.last_wavelet(jmax)) break;
  }
  cout << "  ... done!" << endl;

  cout << "- set up Gramian matrix..." << endl;
  clock_t tstart, tend;
  double time;

  tstart = clock();
  SparseMatrix<double> A;
  setup_stiffness_matrix(problem, Lambda, A);
  tend = clock();

  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- Gramian matrix A=" << endl << A << endl;

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
  cout << "- prepare for the computation of pointwise errors..." << endl;
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
  
  InfiniteVector<double,Index> ucoeffs;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    ucoeffs.set_coefficient(*it, x[i]);

  cout << "- compute L_infty error on a subgrid..." << endl;
  SampledMapping<2> rf(grid, *u);
  SampledMapping<2> rf0(basis.evaluate(ucoeffs, resi));
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
#endif
  
#if 0
  // infinite loop to measure the memory consumption
  while (true) {}
#endif

  if (u) delete u;

  return 0;
}

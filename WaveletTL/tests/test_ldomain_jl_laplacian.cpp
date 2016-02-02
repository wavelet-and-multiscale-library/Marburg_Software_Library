#include <iostream>
#include <set>

#include <utils/function.h>
#include <algebra/sparse_matrix.h>
#include <numerics/iteratsolv.h>
#include <numerics/bezier.h>
#include <numerics/corner_singularity.h>
#include <Ldomain/ldomain_jl_basis.h>
#include <Ldomain/ldomain_jl_evaluate.h>
#include <Ldomain/ldomain_jl_expansion.h>

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#include <galerkin/ldomain_jl_laplacian.h>
#include <galerkin/galerkin_utils.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing LDomainJLLaplacian ..." << endl;

  typedef LDomainJLBasis Basis;
  typedef Basis::Index Index;

  Basis basis;
  
  const int solution = 3;
  Function<2> *uexact = 0, *f = 0;
  switch(solution) {
  case 1:
    uexact = new PolySolution();
    f      = new PolyRHS();
    break;
  case 2:
    uexact = new EigenSolution();
    f      = new EigenRHS();
    break;
  case 3:
    uexact = new CornerSingularity   (Point<2>(0,0), 0.5, 1.5);
    f      = new CornerSingularityRHS(Point<2>(0,0), 0.5, 1.5);
    break;
  default:
    break;
  }

#if 0
  // temporary hack, choose f=1 just to see the inner products
  delete f;
  f = new ConstantFunction<2>(Vector<double>(1, "1"));
#endif

  const int jmax = 2;
  
  typedef LDomainJLLaplacian Problem;
  Problem problem(basis, InfiniteVector<double,Index>());

  InfiniteVector<double,Index> fcoeffs;
  expand(f, basis, true, jmax, fcoeffs);
  problem.set_rhs(fcoeffs);

//   cout << "- integrals of f against the primal wavelets:" << endl
//        << fcoeffs << endl;

  set<Index> Lambda;
  for (Index lambda = basis.first_generator(basis.j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == basis.last_wavelet(jmax)) break;
  }
  
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
#if 0
  A.matlab_output("LdomainJL_laplacian", "A", 1);
#endif


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
  cout << "- extract point values of the solution:" << endl;
  InfiniteVector<double,Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);
  u.scale(&problem, -1);
//   Array1D<SampledMapping<2> > s(evaluate(problem.basis(), u, 1<<5));
  const int N = 200;
  Array1D<SampledMapping<2> > s(evaluate(problem.basis(), u, N));
//   std::ofstream u_Lambda_stream("u_lambda.m");
//   matlab_output(u_Lambda_stream, s);
//   u_Lambda_stream.close();
//   cout << "  ... done, see file 'u_lambda.m'" << endl;

//   for (int i = 0; i <= 2; i++) {
//     s[i].add(-1.0, SampledMapping<2>(s[i], *uexact));
//     cout << "  ... pointwise error on patch " << i << ": " << row_sum_norm(s[i].values()) << endl;
//   }

  double L2error = 0;
  for (int i = 0; i <= 2; i++) {
    s[i].add(-1.0, SampledMapping<2>(s[i], *uexact));
    const double frob = frobenius_norm(s[i].values());
    L2error += frob*frob;
  }
  L2error = sqrt(L2error/(N*N));
  cout << "  ... L_2 error " << L2error << endl;
#endif

  if (uexact) delete uexact;
  if (f) delete f;

  return 0;
}

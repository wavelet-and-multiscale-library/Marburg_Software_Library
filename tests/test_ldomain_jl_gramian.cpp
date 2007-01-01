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
#include <galerkin/ldomain_jl_gramian.h>
#include <galerkin/galerkin_utils.h>

#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing LDomainJLGramian ..." << endl;

  typedef LDomainJLBasis Basis;
  typedef Basis::Index Index;

  Basis basis;
  
  const int solution = 3;
  Function<2> *uexact = 0, *f = 0;
  switch(solution) {
  case 1:
    uexact = new PolySolution();
    f      = new PolySolution();
    break;
  case 2:
    uexact = new EigenSolution();
    f      = new EigenSolution();
    break;
  case 3:
    uexact = new CubicHermiteInterpolant2D_td(2, 0, 0, -2, 1);
    f      = new CubicHermiteInterpolant2D_td(2, 0, 0, -2, 1);
//     uexact = new CubicHermiteInterpolant2D_td(1, 0, 0, -1, 1);
//     f      = new CubicHermiteInterpolant2D_td(1, 0, 0, -1, 1);
//     uexact = new CubicHermiteInterpolant2D_td(1, 0, 0, -1, -1);
//     f      = new CubicHermiteInterpolant2D_td(1, 0, 0, -1, -1);
//     uexact = new CubicHermiteInterpolant2D_td(1, 1, 1, -1, -2);
//     f      = new CubicHermiteInterpolant2D_td(1, 1, 1, -1, -2);
//     uexact = new CubicHermiteInterpolant2D_td(1, 1, 1, 2, 0);
//     f      = new CubicHermiteInterpolant2D_td(1, 1, 1, 2, 0);
    break;
  case 4:
    {
      Point<2> origin(0., 0.);
      uexact = new CornerSingularity   (origin, 0.5, 1.5);
      f      = new CornerSingularity   (origin, 0.5, 1.5);
      break;
    }
    break;
  default:
    break;
  }

#if 1
  {
    cout << "- plot exact solution into a file..." << endl;
    const bool matlab = false; // false -> Octave output
    Grid<2> grid(Point<2>(-1.0, 0.0), Point<2>(0.0, 1.0), 1<<5, 1<<5); // only patch 0
    SampledMapping<2> u_sm(grid, *uexact);
    std::ofstream u_fs("uexact.m");
    
    if (matlab)
      u_sm.matlab_output(u_fs);
    else
      u_sm.octave_output(u_fs);
    
    if (matlab)
      u_fs << "surf(x,y,z)" << endl;
    else
      u_fs << "mesh(x,y,z)" << endl;
    
    u_fs.close();
    cout << "  ... done, see uexact.m" << endl;
  }
#endif

#if 0
  // temporary hack, choose f=1 just to see the inner products
  delete f;
  f = new ConstantFunction<2>(Vector<double>(1, "1"));
#endif

  const int jmax = 4;
  
  typedef LDomainJLGramian Problem;
  Problem problem(basis, InfiniteVector<double,Index>());

  cout << "- expand right-hand side..." << endl;
  InfiniteVector<double,Index> fcoeffs;
  expand(f, basis, true, jmax, fcoeffs);
  problem.set_rhs(fcoeffs);
  cout << "  ... done!" << endl;

//   cout << "- integrals of f against the primal wavelets:" << endl
//        << fcoeffs << endl;

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
#if 0
  A.matlab_output("LdomainJL_gramian", "G", 1);
#endif

#if 0
  cout << "- validate entries of the (unpreconditioned) Gramian matrix:" << endl;

  for (set<Index>::const_iterator itlambda = Lambda.begin(); itlambda != Lambda.end(); ++itlambda) {
    cout << "* checking row " << *itlambda << "..." << endl;
    for (set<Index>::const_iterator itmu = Lambda.begin(); itmu != Lambda.end(); ++itmu) {
      const int N = 2000;
      const double h = 1./N;

      double s = 0;

      FixedArray1D<Array1D<double>,2> lambdavalues, muvalues;
      FixedArray1D<Array1D<double>,2> knots;
      knots[0].resize(N+1);
      knots[1].resize(N+1);

      const Index& lambda = *itlambda;
      const Index& mu = *itmu;
      
      // patch Omega_0 = [-1,0]x[0,1]
      if (lambda.k()[0] <= 0 && lambda.k()[1] >= 0 && mu.k()[0] <= 0 && mu.k()[1] >= 0) {
	for (int k0 = 0; k0 < N; k0++)
	  knots[0][k0] = -1.0+k0*h;
	evaluate(0, lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0], knots[0], lambdavalues[0]);
	evaluate(0, mu.j(), mu.e()[0], mu.c()[0], mu.k()[0], knots[0], muvalues[0]);
	for (int k1 = 0; k1 < N; k1++)
	  knots[1][k1] = k1*h;
	evaluate(0, lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1], knots[1], lambdavalues[1]);
	evaluate(0, mu.j(), mu.e()[1], mu.c()[1], mu.k()[1], knots[1], muvalues[1]);
	for (int k0 = 0; k0 < N; k0++) {
	  for (int k1 = 1; k1 < N; k1++) {
	    s += lambdavalues[0][k0] * lambdavalues[1][k1]
	      * muvalues[0][k0] * muvalues[1][k1];
	  }
	}
      }
      // patch Omega_1 = [-1,0]x[-1,0]
      if (lambda.k()[0] <= 0 && lambda.k()[1] <= 0 && mu.k()[0] <= 0 && mu.k()[1] <= 0) {
	for (int k0 = 0; k0 < N; k0++)
	  knots[0][k0] = -1.0+k0*h;
	evaluate(0, lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0], knots[0], lambdavalues[0]);
	evaluate(0, mu.j(), mu.e()[0], mu.c()[0], mu.k()[0], knots[0], muvalues[0]);
	for (int k1 = 0; k1 < N; k1++)
	  knots[1][k1] = -1.0+k1*h;
	evaluate(0, lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1], knots[1], lambdavalues[1]);
	evaluate(0, mu.j(), mu.e()[1], mu.c()[1], mu.k()[1], knots[1], muvalues[1]);
	for (int k0 = 0; k0 < N; k0++) {
	  for (int k1 = 1; k1 < N; k1++) {
	    s += lambdavalues[0][k0] * lambdavalues[1][k1]
	      * muvalues[0][k0] * muvalues[1][k1];
	  }
	}
      }
      // patch Omega_2 = [0,1]x[-1,0]
      if (lambda.k()[0] >= 0 && lambda.k()[1] <= 0 && mu.k()[0] >= 0 && mu.k()[1] <= 0) {
	for (int k0 = 0; k0 < N; k0++)
	  knots[0][k0] = k0*h;
	evaluate(0, lambda.j(), lambda.e()[0], lambda.c()[0], lambda.k()[0], knots[0], lambdavalues[0]);
	evaluate(0, mu.j(), mu.e()[0], mu.c()[0], mu.k()[0], knots[0], muvalues[0]);
	for (int k1 = 0; k1 < N; k1++)
	  knots[1][k1] = -1.0+k1*h;
	evaluate(0, lambda.j(), lambda.e()[1], lambda.c()[1], lambda.k()[1], knots[1], lambdavalues[1]);
	evaluate(0, mu.j(), mu.e()[1], mu.c()[1], mu.k()[1], knots[1], muvalues[1]);
	for (int k0 = 0; k0 < N; k0++) {
	  for (int k1 = 1; k1 < N; k1++) {
	    s += lambdavalues[0][k0] * lambdavalues[1][k1]
	      * muvalues[0][k0] * muvalues[1][k1];
	  }
	}
      }
      
      const double alambdamu = problem.a(*itlambda,*itmu);
      const double alambdamuapprox = s*h*h;
      const double dev = fabs(alambdamu-alambdamuapprox);
      if (dev > 1e-6)
	cout << "lambda=" << *itlambda
	     << ", mu=" << *itmu << ": "
	     << "a(.,.)=" << alambdamu
	     << ", approx.=" << alambdamuapprox
	     << ", dev.=" << dev << endl;
    }
  }
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
  
#if 0
  {
    cout << "- plot point values of the solution:" << endl;
    InfiniteVector<double,Index> u;
    unsigned int i = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
      u.set_coefficient(*it, x[i]);
    u.scale(&problem, -1);
    Array1D<SampledMapping<2> > s(evaluate(problem.basis(), u, 1<<5));
    std::ofstream u_Lambda_stream("u_lambda.m");
    octave_output(u_Lambda_stream, s);
    u_Lambda_stream.close();
    cout << "  ... done, see file 'u_lambda.m'" << endl;
    
    for (int i = 0; i <= 2; i++) {
      s[i].add(-1.0, SampledMapping<2>(s[i], *uexact));
      cout << "  pointwise error on patch " << i << ": " << row_sum_norm(s[i].values()) << endl;
    }
  }
#endif

#if 1
  {
    cout << "- compute error:" << endl;
    InfiniteVector<double,Index> u;
    unsigned int i = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
      u.set_coefficient(*it, x[i]);
    u.scale(&problem, -1);
    const int N = 64;
    Array1D<SampledMapping<2> > s(evaluate(problem.basis(), u, N));
    std::ofstream u_Lambda_stream("u_lambda.m");
    octave_output(u_Lambda_stream, s);
    u_Lambda_stream.close();
    cout << "  ... plot of the approximation: see u_lambda.m" << endl;
    
    double L2error = 0;
    for (int ii = 0; ii <= 2; ii++) {
      s[ii].add(-1.0, SampledMapping<2>(s[ii], *uexact));
      const double frob = frobenius_norm(s[ii].values());
      L2error += frob*frob;
    }
    L2error = sqrt(L2error/(N*N));
    cout << "  ... L_2 error " << L2error << endl;

    cout << "- plot error into a file:" << endl;
    std::ofstream u_Lambda_error_stream("u_lambda_error.m");
    octave_output(u_Lambda_error_stream, s);
    u_Lambda_error_stream.close();
    cout << "  ... done, see file 'u_lambda_error.m'" << endl;
  }
#endif

  if (uexact) delete uexact;
  if (f) delete f;

  return 0;
}

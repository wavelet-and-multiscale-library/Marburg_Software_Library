#include <iostream>
#include <set>

#include <utils/function.h>
#include <algebra/sparse_matrix.h>
#include <interval/spline_basis.h>
#include <numerics/iteratsolv.h>
#include <numerics/corner_singularity.h>
#define _WAVELETTL_LDOMAINBASIS_VERBOSITY 0
#include <Ldomain/ldomain_basis.h>
#include <Ldomain/ldomain_evaluate.h>
#include <Ldomain/ldomain_expansion.h>

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#include <galerkin/galerkin_utils.h>
#include <galerkin/ldomain_helmholtz_equation.h>
#include "ldomain_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing LDomainHelmholtzEquation ..." << endl;

  const int d  = 2;
  const int dT = 2;

  typedef SplineBasis<d,dT,DS_construction> Basis1D;
  Basis1D basis1d("bio5-energy",0,0,0,0);
  
  typedef LDomainBasis<Basis1D> Basis;
  typedef Basis::Index Index;
  Basis basis(basis1d);
  
  const int solution = 1;
  Function<2> *uexact = 0, *f = 0;
  switch(solution) {
  case 1:
    uexact = new PolySolution();
    f = new PolySum(); // alpha=1
//     f = new PolyRHS(); // alpha=0
    break;
  default:
    break;
  }

//   const int jmax = basis.j0()+1;
  const int jmax = 5;
  
  typedef LDomainHelmholtzEquation<d,dT> Problem;
  Problem problem(basis,
		  "DS_B_2_2_5_G",
		  "DS_B_2_2_5_A",
		  5,
		  1.0,
		  InfiniteVector<double, Index>());

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
  A.matlab_output("Ldomain_laplacian", "A", 1);
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
  {
    cout << "- plot point values of the solution:" << endl;
    InfiniteVector<double,Index> u;
    unsigned int i = 0;
    for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
      u.set_coefficient(*it, x[i]);
    u.scale(&problem, -1);
    const int resi = 5;
    const int N = 1<<resi;
    Array1D<SampledMapping<2> > s(basis.evaluate(u, resi));
    std::ofstream u_Lambda_stream("u_lambda.m");
//     octave_output(u_Lambda_stream, s);
    matlab_output(u_Lambda_stream, s);
    u_Lambda_stream.close();
    cout << "  ... done, see file 'u_lambda.m'" << endl;
 
    double L2error = 0;
    for (int i = 0; i <= 2; i++) {
      s[i].add(-1.0, SampledMapping<2>(s[i], *uexact));
      cout << "  pointwise error on patch " << i << ": " << row_sum_norm(s[i].values()) << endl;
      const double frob = frobenius_norm(s[i].values());
      L2error += frob*frob;
    }
    L2error = sqrt(L2error/(N*N));
    cout << "  ... L_2 error " << L2error << endl;
  }
#endif

  if (uexact) delete uexact;
  if (f) delete f;

  return 0;
}

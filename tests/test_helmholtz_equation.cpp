#include <iostream>
#include <set>
#include <interval/spline_basis.h>
#include <interval/spline_expansion.h>
#include <galerkin/helmholtz_equation.h>
#include <utils/function.h>
#include <numerics/iteratsolv.h>

#include "helmholtz_1d_solutions.h"

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

int main()
{
  cout << "Testing HelmholtzEquation ..." << endl;

  const int d  = 3;
  const int dT = 3;

  typedef SplineBasis<d,dT,P_construction,1,1,0,0> Basis; // PBasis, complementary b.c.'s
  Basis basis;  

  typedef Basis::Index Index;


  const unsigned int solution = 1;
  double kink = 0; // for Solution3;

  Function<1> *uexact = 0, *f = 0;
  switch(solution) {
  case 1:
    uexact = new Solution1();
    f = new RHS1();
    break;
  case 2:
    uexact = new Solution2();
    f = new RHS2();
    break;
  case 3:
    kink = 5./7.;
    uexact = new Solution3(kink);
    f = new RHS3_part(kink);
    break;
  default:
    break;
  }

  int jmax = 9;
  basis.set_jmax(jmax);


  typedef HelmholtzEquation1D<d,dT> Problem;
  Problem helmholtz(basis, 1.0, InfiniteVector<double,Index>());

  InfiniteVector<double,Index> fcoeffs;
  expand(f, basis, true, jmax, fcoeffs);
  helmholtz.set_rhs(fcoeffs);

  InfiniteVector<double, Index> coeffs;
  helmholtz.RHS(1e-8, coeffs);

  jmax = helmholtz.basis().j0()+1;

  set<Index> Lambda;
  for (Index lambda = helmholtz.basis().first_generator(helmholtz.basis().j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == helmholtz.basis().last_wavelet(jmax-1)) break;
  }
  
//   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//     cout << *it << endl;

  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();
  
  SparseMatrix<double> A;
  setup_stiffness_matrix(helmholtz, Lambda, A);

  cout << "- stiffness matrix: " << endl << A << endl;
#if 0
  std::ofstream ofs("stiff_out.m");
  ofs << "M=";
  print_matrix(A,ofs);
  ofs.close();
#endif

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;

  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(helmholtz, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- right hand side: " << b << endl;

  Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
  unsigned int iterations;
  CG(A, b, x, 1e-15, 100, iterations);
  
//   cout << "- solution coefficients: " << x;
  cout << " with residual (infinity) norm ";
  A.apply(x, err);
  err -= b;
  cout << linfty_norm(err) << endl;

  // apply D^{-1} to obtain L_2 wavelet coefficients
  InfiniteVector<double,Index> ulambda;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    ulambda.set_coefficient(*it, x[i]);
  ulambda.scale(&helmholtz, -1);

  // insert coefficients into a dense vector
  Vector<double> wcoeffs(helmholtz.basis().Deltasize(jmax));
  for (InfiniteVector<double,Index>::const_iterator it(ulambda.begin()),
	 itend(ulambda.end()); it != itend; ++it) {
    // determine number of the wavelet
    typedef Vector<double>::size_type size_type;
    size_type number = 0;
    if (it.index().e() == 0) {
      number = it.index().k()-helmholtz.basis().DeltaLmin();
    } else {
      number = helmholtz.basis().Deltasize(it.index().j())+it.index().k()-helmholtz.basis().Nablamin();
    }
    wcoeffs[number] = *it;
  }
  
  // switch to generator representation
  Vector<double> gcoeffs(wcoeffs.size(), false);
  if (jmax == helmholtz.basis().j0())
    gcoeffs = wcoeffs;
  else
    helmholtz.basis().apply_Tj(jmax-1, wcoeffs, gcoeffs);
  
  // evaluate Galerkin solution
  const unsigned int N = 100;
  const double h = 1./N;
  Vector<double> ulambda_values(N+1);
  for (unsigned int i = 0; i <= N; i++) {
    const double x = i*h;
    SchoenbergIntervalBSpline_td<d> sbs(jmax,0);
    for (unsigned int k = 0; k < gcoeffs.size(); k++) {
      sbs.set_k(helmholtz.basis().DeltaLmin()+k);
      ulambda_values[i] += gcoeffs[k] * sbs.value(Point<1>(x));
    }
  }

  // evaluate exact solution
  Vector<double> uexact_values(N+1);
  for (unsigned int i = 0; i <= N; i++) {
    const double x = i*h;
    uexact_values[i] = uexact->value(Point<1>(x));
  }
  
  const double Linfty_error = linfty_norm(ulambda_values-uexact_values);
  cout << "  L_infinity error on a subgrid: " << Linfty_error << endl;


  cout << "The stiffness matrix again: A=" << endl << A << endl;
  //InfiniteVector<double,Index> w;
  Vector<double> w(basis.degrees_of_freedom());
  Index lambda = helmholtz.basis().first_generator(helmholtz.basis().j0());
  lambda = helmholtz.basis().first_wavelet(helmholtz.basis().j0());
  ++lambda;
  cout << "lambda=" << lambda << ", call add_level():" << endl;
  helmholtz.add_level(lambda, w, helmholtz.basis().j0(), 1.0, 1000);
  cout << w << endl;
  

  if (f) delete f;
  if (uexact) delete uexact;

  return 0;
}

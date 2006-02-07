#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <galerkin/sturm_equation.h>

using namespace std;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    switch(N) {
    case 1:
      return 1;
      break;
    default:
      return 0;
      break;
    }
  }
  double p_prime(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double q(const double t) const {
    switch(N) {
    case 1:
      return 0;
      break;
    default:
      return 0;
      break;
    }
  }
  double g(const double t) const {
    switch(N) {
    case 1:
      return 2;
      break;
    default:
      return 0;
      break;
    }
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};


int main()
{
  cout << "Testing wavelet-Galerkin solution of a Sturm b.v.p. ..." << endl;

  TestProblem<1> T;

  const int d  = 2;
  const int dT = 2; // be sure to use a continuous dual here, otherwise the RHS test will fail
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;

  SturmEquation<Basis> eq(T);

  InfiniteVector<double, Index> coeffs;

#if 0
  coeffs[first_generator(&eq.basis(), eq.basis().j0())] = 1.0;
  coeffs[last_generator(&eq.basis(), eq.basis().j0())] = 2.0;
  coeffs[first_wavelet(&eq.basis(), eq.basis().j0())] = 3.0;
  coeffs[last_wavelet(&eq.basis(), eq.basis().j0())] = 4.0;
  coeffs[first_wavelet(&eq.basis(), eq.basis().j0()+1)] = 5.0;
  cout << "- a coefficient set:" << endl
       << coeffs << endl;
  eq.rescale(coeffs, -1);
  cout << "- after rescaling with D^{-1}:" << endl
       << coeffs << endl;
#endif

  eq.RHS(1e-2, coeffs);

#if 0
  cout << "- approximate coefficient set of the right-hand side:" << endl
       << coeffs << endl;
  cout << "- check expansion of the right-hand side in the dual basis:" << endl;
  eq.rescale(coeffs, 1);
  evaluate(eq.basis(), coeffs, false, 8).matlab_output(cout);
  eq.rescale(coeffs, -1);
#endif  

  set<Index> Lambda;
  for (Index lambda = first_generator(&eq.basis(), eq.basis().j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == last_wavelet(&eq.basis(), eq.basis().j0())) break;
  }

//   cout << "- set up stiffness matrix with respect to the index set Lambda=" << endl;
//   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it)
//     cout << *it << endl;

#if 1
  cout << "- set up (preconditioned) stiffness matrix..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();
  SparseMatrix<double> A;
  setup_stiffness_matrix(eq, Lambda, A);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
//   cout << "- (preconditioned) stiffness matrix A=" << endl << A << endl;

  cout << "- set up right-hand side..." << endl;
  tstart = clock();
  Vector<double> b;
  setup_righthand_side(eq, Lambda, b);
  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;
  cout << "- right hand side: " << b << endl;

  Vector<double> x(Lambda.size()), err(Lambda.size()); x = 0;
  unsigned int iterations;
  CG(A, b, x, 1e-8, 100, iterations);
  
  cout << "- solution coefficients: " << x;
  cout << " with residual (infinity) norm ";
  A.apply(x, err);
  err -= b;
  cout << linfty_norm(err) << endl;
  
  cout << "- point values of the solution:" << endl;
  InfiniteVector<double,Index> u;
  unsigned int i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, x[i]);
  
  u.scale(&eq, -1);
  SampledMapping<1> s(evaluate(eq.basis(), u, true, 7));
  s.matlab_output(cout);
#endif

#if 0
  cout << "- estimate for ||D^{-1}AD^{-1}||: " << eq.norm_A() << endl;
  cout << "- estimate for ||(D^{-1}AD^{-1})^{-1}||: " << eq.norm_Ainv() << endl;
#endif

#if 0
  cout << "- checking add_column and the compression strategy:" << endl;
  InfiniteVector<double,Index> v, w;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it) {
    v.clear();
    w.clear();
    eq.add_column(1.0, *it, 99, w);
    cout << "lambda=" << *it << endl;
//     cout << "* result of add_column:" << endl << w << endl;
    const int jmax = 8;
    cout << "* checking correctness up to level " << jmax << "..." << endl;
    for (Index nu = first_generator(&eq.basis(), eq.basis().j0());; ++nu) {
      v.set_coefficient(nu, eq.a(nu, *it)/(eq.D(nu)*eq.D(*it)));
      if (nu == last_wavelet(&eq.basis(), jmax)) break;
    }
//     cout << "  ... done, error:" << endl;
//     InfiniteVector<double,Index> error = v - w;
//     error.compress();
//     cout << error << endl;

    cout << "  ... done, absolute error " << linfty_norm(v-w) << endl;

    break;
  }

  cout << "- check add_column for a second time to test the cache:" << endl;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it) {
    v.clear();
    w.clear();
    eq.add_column(1.0, *it, 99, w);
    cout << "lambda=" << *it << endl;
//     cout << "* result of add_column:" << endl << w << endl;
    const int jmax = 8;
    cout << "* checking correctness up to level " << jmax << "..." << endl;
    for (Index nu = first_generator(&eq.basis(), eq.basis().j0());; ++nu) {
      v.set_coefficient(nu, eq.a(nu, *it)/(eq.D(nu)*eq.D(*it)));
      if (nu == last_wavelet(&eq.basis(), jmax)) break;
    }
    cout << "  ... done, absolute error " << linfty_norm(v-w) << endl;

    break;
  }
#endif

  InfiniteVector<double,Index> coeffsp;
  eq.left_preconditioner()->apply_preconditioner(coeffs,coeffsp);

  return 0;
}

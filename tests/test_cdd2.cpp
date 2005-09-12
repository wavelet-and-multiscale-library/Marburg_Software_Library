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

#include <adaptive/cdd2.h>

using namespace std;
using namespace WaveletTL;

using MathTL::simpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
 */
template <unsigned int N>
class TestProblem
  : public simpleSturmBVP
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
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD2_SOLVE ..." << endl;

  TestProblem<1> T;

  const int d  = 2;
  const int dT = 4;
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;

  SturmEquation<Basis> problem(T);

  InfiniteVector<double, Index> F_eta;
  problem.RHS(1e-6, F_eta);
  const double nu = problem.norm_Ainv() * l2_norm(F_eta);

  InfiniteVector<double, Index> u_epsilon;
  CDD2_SOLVE(problem, nu, 1e-1, u_epsilon);
  
  return 0;
}

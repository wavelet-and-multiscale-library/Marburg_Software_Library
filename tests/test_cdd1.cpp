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
#include <galerkin/cached_problem.h>

#define _WAVELETTL_CDD1_VERBOSITY 1
#include <adaptive/cdd1.h>

using namespace std;
using namespace WaveletTL;

using MathTL::SimpleSturmBVP;
using MathTL::CG;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  1: y(t)=x*(1-x), -y''(t)=2
  2: y(t)=exp(-100*(x-0.5)^2), -y''(t)= (200-(200x-100)^2)*exp(-100*(x-0.5)^2)
  3: 1D example from [BBCCDDU]
 */
template <unsigned int N>
class TestProblem
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    switch(N) {
    case 1:
    case 2:
    case 3:
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
    case 2:
    case 3:
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
    case 2:
    case 3:
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
    case 2:
      return (200.0-(200.0*t-100.0)*(200.0*t-100.0))*exp(-100.0*(t-0.5)*(t-0.5));
      break;
    case 3:
      return -100*exp(5*t)*(1-(exp(5*t)-1)/(exp(5.)-1))/(exp(5.)-1)+200*exp(10*t)/((exp(5.)-1)*(exp(5.)-1))+100*(exp(5*t)-1)*exp(5*t)/((exp(5.)-1)*(exp(5.)-1));
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
  cout << "Testing adaptive wavelet-Galerkin solution of a Sturm b.v.p. with CDD1_SOLVE ..." << endl;

  TestProblem<3> T;

  const int d  = 3;
  const int dT = 3;
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;

  SturmEquation<Basis> problem(T);
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem);
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem,  2.34801,  13.3113 ); // d=2, dT=2 (diag-precond.)
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem, 12.2509 ,   6.41001); // d=2, dT=2 (2^j-precond.)
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem,  3.04627,  52.2117 ); // d=2, dT=4 (diag-precond.)
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem, 15.6751 ,  25.9767 ); // d=2, dT=4 (2^j-precond.)
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem,  2.12054, 209.511  ); // d=3, dT=3 (diag-precond.)
  CachedProblem<SturmEquation<Basis> > cproblem(&problem,  6.7774 ,  45.576  ); // d=3, dT=3 (2^j-precond.)
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem,  2.92919, 581.862  ); // d=3, dT=5 (diag-precond.)
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem, 13.5193 , 103.532  ); // d=3, dT=5 (2^j-precond.)
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem, 38.4615 , 251.849  ); // d=3, dT=7
//   CachedProblem<SturmEquation<Basis> > cproblem(&problem,  5.8466 , 532.085  ); // d=4, dT=4
  cout << "* estimate for normA: " << cproblem.norm_A() << endl;
  cout << "* estimate for normAinv: " << cproblem.norm_Ainv() << endl;

  InfiniteVector<double, Index> u_epsilon;
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon, 8);
  CDD1_SOLVE(cproblem, 1e-2, u_epsilon, 10);
//   CDD1_SOLVE(cproblem, 1e-4, u_epsilon);
  
  return 0;
}

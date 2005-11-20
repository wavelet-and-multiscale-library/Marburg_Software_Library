#include <iostream>
#include <fstream>
#include <set>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>
#include <parabolic/lin_par_equation.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

/*
  different test problems with homogeneous Dirichlet b.c.'s
  0: y=0, y''(t)=0
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
    case 0:
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
    case 0:
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
    case 0:
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
    case 0:
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
  cout << "Testing linear parabolic equations..." << endl;

  TestProblem<0> T;
  
  const int d  = 2;
  const int dT = 2;
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;
  typedef SturmEquation<Basis> EllipticEquation;
  typedef InfiniteVector<double,Index> V;

  EllipticEquation elliptic(T);
//   CachedProblem<EllipticEquation> celliptic(&elliptic);
  CachedProblem<EllipticEquation> celliptic(&elliptic, 12.2508, 6.41001); // d=2, dT=2

  V u0;
  u0.set_coefficient(Index(3,0,4, &celliptic.basis()), 1.0); // an inner hat fct.

  LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0);
  
  cout << "* testing ROS2:" << endl;
  ROWMethod<V> ros2(WMethod<V>::ROS2);
//   scheme = &ros2;
//   olderr = 0;
//   for (int expo = 0; expo <= 6; expo++) {
//     temp = problem.u0;
//     int N = 1<<expo;
//     double h = 1.0/N;
//     double err_est = 0;
//     for (int i = 1; i <= N; i++) {
//       scheme->increment(&problem, i*h, temp, h, result, error_estimate);
//       err_est = std::max(err_est, l2_norm(error_estimate));
//       temp = result;
//     }
//     err = fabs(result[0] - problem.exact_solution(1.0));
//     if (expo > 0) {
//       cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
// 	   << ", max. of local error estimate " << err_est << endl;
//     }
//     olderr = err;
//   }


  return 0;
}

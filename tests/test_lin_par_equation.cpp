#include <iostream>
#include <fstream>
#include <set>
#include <list>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <interval/ds_expansion.h>
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

class Hat
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return std::max(0.0,0.5-abs(p[0]-0.5));
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class Haar
  : public Function<1>
{
  public:
  inline double value(const Point<1>& p,
		      const unsigned int component = 0) const
  {
    return abs(p[0]-0.5)<=0.25 ? 1.0 : 0.0;
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing linear parabolic equations..." << endl;

  TestProblem<0> testproblem;
  
  const int d  = 2;
  const int dT = 2;
  typedef DSBasis<d,dT> Basis;
  typedef Basis::Index Index;
  typedef SturmEquation<Basis> EllipticEquation;
  typedef InfiniteVector<double,Index> V;

  EllipticEquation elliptic(testproblem);
//   CachedProblem<EllipticEquation> celliptic(&elliptic);
  CachedProblem<EllipticEquation> celliptic(&elliptic, 12.2508, 6.41001); // d=2, dT=2
//   CachedProblem<SturmEquation<Basis> > celliptic(&elliptic, 6.73618, 45.5762); // d=3, dT=3

//   Hat hat;
  Haar hat;
  V u0;
  expand(&hat, celliptic.basis(), false, celliptic.basis().j0()+4, u0);
  u0.compress(1e-14);
  cout << "* expansion coefficients of hat function u0=" << endl << u0;

  LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0);
  
#if 0
  cout << "* testing ROS2:" << endl;
  ROWMethod<V> ros2(WMethod<V>::ROS2);
  OneStepScheme<V>* scheme = &ros2;
  V temp, result, error_estimate;
  double err, olderr = 0;

  olderr = 0;
//   for (int expo = 0; expo <= 6; expo++) {
  for (int expo = 2; expo <= 2; expo++) {
    temp = parabolic.u0;
    
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      cout << "---------------- before increment() -----------------------" << endl;
      scheme->increment(&parabolic, i*h, temp, h, result, error_estimate, 1e-1);
      cout << "---------------- after increment() -----------------------" << endl;
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
//     err = fabs(result[0] - problem.exact_solution(1.0));
//     if (expo > 0) {
    cout << "h=" << h << ", max. of local error estimate " << err_est << endl;
//     }
//     olderr = err;
  }
#endif

#if 1
  const double T = 1.0;
  const double q = 10.0;
//   const double TOL = 1e-2;
  const double TOL = 1e-1;
  const double tau_max = 1.0;

  IVPSolution<V> result_adaptive;
  
  cout << "* testing ROS2 (adaptive)..." << endl;
  ROWMethod<V> ros2_adaptive(WMethod<V>::ROS2);
  solve_IVP(&parabolic, &ros2_adaptive, T,
	    TOL, q, tau_max, result_adaptive);

  const int resolution = 10;
  Grid<1> gx(0.0, 1.0, 1<<resolution);
  Array1D<double> points(result_adaptive.t.size());
  std::copy(result_adaptive.t.begin(), result_adaptive.t.end(), points.begin());
  Grid<1> gt(points);
  Grid<2> grid(gt, gx);
  Matrix<double> values((1<<resolution)+1, result_adaptive.t.size());
  
  unsigned int i(0);
  for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
       ui != result_adaptive.u.end(); ++ui, ++i) {
    SampledMapping<1> temp(evaluate(celliptic.basis(), *ui, true, resolution));
    for (unsigned int k(0); k < values.row_dimension(); k++)
      values.set_entry(k, i, temp.values()[k]);
  }
  SampledMapping<2> result(grid, values);
 
  std::ofstream resultstream("lin_par_result.m");
  result.matlab_output(resultstream);
  resultstream.close();
 
#endif

  return 0;
}

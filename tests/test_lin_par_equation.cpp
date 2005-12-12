#include <iostream>
#include <fstream>
#include <set>
#include <list>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>
#define _MATHTL_ONESTEPSCHEME_VERBOSITY 1
#include <numerics/w_method.h>
#include <numerics/row_method.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <interval/ds_expansion.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>

#define _WAVELETTL_CDD1_VERBOSITY 0
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

class Hat : public Function<1>
{
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return std::max(0.0,0.5-abs(p[0]-0.5));
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class Haar : public Function<1>
{
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return abs(p[0]-0.5)<=0.25 ? 1.0 : 0.0;
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class constant_f : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return M_PI*M_PI*sin(M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class rhs_4 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return -5*p[0]+9*p[0]*p[0]+2*p[0]*p[0]*p[0]+6*get_time()-12*get_time()*p[0]; // Maple...
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class rhs_5 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return p[0]*(1-p[0])+2*get_time();
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class exact_solution_3 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return (1-exp(-M_PI*M_PI*get_time()))*sin(M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class exact_solution_4 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return get_time()*p[0]*(1-p[0])*(1-p[0])*(1-p[0]) + (1-get_time())*p[0]*p[0]*p[0]*(1-p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class exact_solution_5 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return get_time()*p[0]*(1-p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
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

  EllipticEquation elliptic(testproblem, false); // do not precompute the dummy rhs

//   CachedProblem<EllipticEquation> celliptic(&elliptic);
  CachedProblem<EllipticEquation> celliptic(&elliptic, 12.2508, 6.41001); // d=2, dT=2
//   CachedProblem<SturmEquation<Basis> > celliptic(&elliptic, 6.73618, 45.5762); // d=3, dT=3

  const int jmax = 10;

  // handle different test cases:
  // 1: u0 = hat function, f(t)=0
  // 2: u0 = haar function, f(t)=0
  // 3: u0 = 0, f(t,x) = pi^2*sin(pi*x), u(t,x)=(1-exp(-pi^2*t))*sin(pi*x)
  // 4: u(t,x)=t*x*(1-x)^3+(1-t)*x^3*(1-x), f(t,x)=u_t(t,x)-u_{xx}(t,x)
  // 5: u(t,x)=t*x*(1-x), f(t,x) = x*(1-x)+2*t

#define _TESTCASE 3

  //
  //
  // setup exact solution (if available)

#if _TESTCASE == 3
  exact_solution_3 uexact;
#endif

#if _TESTCASE == 4
  exact_solution_4 uexact;
#endif

#if _TESTCASE == 5
  exact_solution_5 uexact;
#endif

  //
  //
  // setup initial value u0

  V u0;

#if _TESTCASE == 1
  Hat hat;
  expand(&hat, celliptic.basis(), false, jmax, u0); 
#endif
  
#if _TESTCASE == 2
  Haar haar;
  expand(&haar, celliptic.basis(), false, jmax, u0); 
#endif

#if _TESTCASE == 3
  // do nothing, u0=0
#endif

#if _TESTCASE == 4
  uexact.set_time(0);
  expand(&uexact, celliptic.basis(), false, jmax, u0);
#endif

#if _TESTCASE == 5
  // do nothing, u0=0
#endif

  u0.compress(1e-14);
//   cout << "* expansion coefficients of u0=" << endl << u0;
 
  //
  //
  // setup driving term f

#if _TESTCASE == 1
  LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0, 0, jmax);
#endif

#if _TESTCASE == 2
  LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0, 0, jmax);
#endif

#if _TESTCASE == 3
  V f;
  constant_f cf;
  expand(&cf, celliptic.basis(), false, jmax, f);
  f.compress(1e-14);
//   LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0, f, jmax);
  LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0, &cf, jmax);
#endif

#if _TESTCASE == 4
  rhs_4 f4;
  LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0, &f4, jmax);
#endif
  
#if _TESTCASE == 5
  rhs_5 f5;
  LinearParabolicEquation<CachedProblem<EllipticEquation> > parabolic(&celliptic, u0, &f5, jmax);
#endif
  
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

#if 0
  // einzelner Testlauf, gibt Plot der Iterierten aus

  const double T = 1.0;
  const double q = 10.0;
//   const double TOL = 1e-2;
  const double TOL = 1e-2;
  const double tau_max = 1.0;

  IVPSolution<V> result_adaptive;
  
  cout << "* testing ROS2 (adaptive, single run)..." << endl;
  ROWMethod<V> ros2_adaptive(WMethod<V>::ROS2);
  solve_IVP(&parabolic, &ros2_adaptive, T,
	    TOL, 0, q, tau_max, result_adaptive);

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

#if 1
  // mehrere Testlaeufe mit einem Problem, verschiedene Toleranzen

  const double T = 1.0;
  const double q = 5.0;
  const double tau_max = 1.0;

  std::list<double> numberofsteps;
  std::list<double> errors;

  cout << "* testing ROS2 (adaptive, several tolerances)..." << endl;
  for (int expo = 10; expo <= 12; expo++) { // 2^{-6}=0.015625, 2^{-8}=3.9e-3, 2^{-10}=9.77e-4
    const double TOL = ldexp(1.0, -expo);

    IVPSolution<V> result_adaptive;

    cout << "  TOL=" << TOL << endl;

    // adaptive solution of u'=Au+f
    ROWMethod<V> ros2_adaptive(WMethod<V>::ROS2);
    solve_IVP(&parabolic, &ros2_adaptive, T,
	      TOL, 0, q, tau_max, result_adaptive);

#if _TESTCASE == 3 || _TESTCASE == 4 || _TESTCASE == 5
    // compute maximal ell_2 error of the coefficients
    double errhelp = 0;
    std::list<double>::const_iterator ti(result_adaptive.t.begin());
    for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
	 ui != result_adaptive.u.end(); ++ui, ++ti) {
      V uexact_coeffs;
      uexact.set_time(*ti);
      expand(&uexact, celliptic.basis(), false, jmax, uexact_coeffs);
      cout << "  ell_2 error at t=" << *ti << ": " << l2_norm(*ui - uexact_coeffs) << endl;
      errhelp = std::max(errhelp, l2_norm(*ui - uexact_coeffs));
    }
    errors.push_back(errhelp);
    numberofsteps.push_back(result_adaptive.t.size());
#endif

#if _TESTCASE == 1 || _TESTCASE == 2
    // quick hack: use TOL
    errors.push_back(TOL);
    numberofsteps.push_back(result_adaptive.t.size());
#endif
  }

  std::ofstream resultstream("work_precision.m");

  resultstream << "errors=[";
  for (std::list<double>::const_iterator it = errors.begin();
       it != errors.end(); ++it) {
    resultstream << log10(*it);
    if (it != errors.end())
      resultstream << " ";
  }
  resultstream << "];" << endl;

  resultstream << "N=[";
  for (std::list<double>::const_iterator it = numberofsteps.begin();
       it != numberofsteps.end(); ++it) {
    resultstream << log10(*it);
    if (it != numberofsteps.end())
      resultstream << " ";
  }
  resultstream << "];" << endl;

  resultstream.close();

#endif

  return 0;
}

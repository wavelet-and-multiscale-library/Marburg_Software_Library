#include <cmath>
#include <iostream>
#include <list>
#include <map>
#include <algebra/vector.h>
#include <numerics/ivp.h>

#define _MATHTL_ONESTEPSCHEME_VERBOSITY 1
#include <numerics/one_step_scheme.h>

#include <numerics/runge_kutta.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>

using std::cout;
using std::endl;
using namespace MathTL;

/*
  The Dahlquist standard problem
    u' = lambda*u, u(0) = 1
  
  So f(t,u)=lambda*u, f_t(t,u)=0 and f_u(t,u)=lambda*I.
 */
class Dahlquist
  : public AbstractIVP<Vector<double> >
{
public:
  Dahlquist(const double lambda = 1)
    : lambda_(lambda)
  {
    u0.resize(1); u0[0] = 1;
  }

  void evaluate_f(const double t,
		  const Vector<double>& v,
		  const double tolerance,
		  Vector<double>& result) const
  {
    result[0] = lambda_ * v[0];
  }

  void evaluate_ft(const double t,
		   const Vector<double>& v,
		   const double tolerance,
		   Vector<double>& result) const
  {
    result = 0;
  }

  void solve_ROW_stage_equation(const double t,
				const Vector<double>& v,
				const double alpha,
				const Vector<double>& y,
				const double tolerancs,
				Vector<double>& result) const
  {
    // Jv=lambda*v -> (alpha*I-J)x=(alpha-lambda)x
    result[0] = y[0] / (alpha - lambda_);
  }

  // exact solution
  double exact_solution(const double t) const
  {
    return exp(lambda_*t);
  }

private:
  double lambda_;
};

/*
  The circle u(t)=(cos(t),sin(t)) with
    u' = [0 -1; 1 0]*u(t), u(0) = [1 0]'
  
  So f(t,u)=lambda*u, f_t(t,u)=0 and f_u(t,u)=lambda*I.
 */
class Circle
  : public AbstractIVP<Vector<double> >
{
public:
  Circle()
  {
    u0.resize(2); u0[0] = 1; u0[1] = 0;
  }

  void evaluate_f(const double t,
		  const Vector<double>& v,
		  const double tolerance,
		  Vector<double>& result) const
  {
    result[0] = -v[1];
    result[1] = v[0];
  }

  void evaluate_ft(const double t,
		   const Vector<double>& v,
		   const double tolerance,
		   Vector<double>& result) const
  {
    result = 0;
  }

  void solve_ROW_stage_equation(const double t,
				const Vector<double>& v,
				const double alpha,
				const Vector<double>& y,
				const double tolerancs,
				Vector<double>& result) const
  {
    const double det_alphaIminusA = alpha*alpha + 1.0;
    result[0] = (alpha*y[0] - y[1]) / det_alphaIminusA;
    result[1] = (y[0] + alpha*y[1]) / det_alphaIminusA;
  }

  // exact solution
  void exact_solution(const double t, Vector<double>& y) const
  {
    y.resize(2);
    y[0] = cos(t);
    y[1] = sin(t);
  }

private:
  double lambda_;
};

/*
  Arenstorf orbit

    y_1'' = y_1 + 2*y_2' - mu'*(y_1+mu)/D_1 - mu*(y_1-mu')/D_2
    y_2'' = y_2 - 2*y_1' - mu'*y_2/D_1 - mu*y_2/D_2

  where

    D_1 = ((y_1+mu)^2+y_2^2)^{3/2},
    D_2 = ((y_1-mu')^2+y_2^2)^{3/2},
    mu = 0.012277471, mu' = 1 - mu

  With initial values y_1(0)=0.994, y_1'(0) = y_2(0) = 0 and
  y_2'(0) = -2.00158510637908252240537862224, the solution is
  periodic, x_end = 17.0652165601579625588917206249
 */
class Arenstorf
  : public AbstractIVP<Vector<double> >
{
public:
  Arenstorf()
  {
    u0.resize(4);
    u0[0] = 0.994;
    u0[1] = u0[2] = 0;
    u0[3] = -2.00158510637908252240537862224;
  }

  void evaluate_f(const double t,
		  const Vector<double>& v,
		  const double tolerance,
		  Vector<double>& result) const
  {
    const double mu = 0.012277471;

    result[0] = v[1];
    result[1] = v[0] + 2*v[3]
      - (1-mu)*(v[0]+mu)/pow((v[0]+mu)*(v[0]+mu)+v[2]*v[2], 1.5)
      - mu*(v[0]-1+mu)/pow((v[0]-1+mu)*(v[0]-1+mu)+v[2]*v[2], 1.5);
    result[2] = v[3];
    result[3] = v[2] - 2*v[1]
      - (1-mu)*v[2]/pow((v[0]+mu)*(v[0]+mu)+v[2]*v[2], 1.5)
      - mu*v[2]/pow((v[0]-1+mu)*(v[0]-1+mu)+v[2]*v[2], 1.5);
  }

  void evaluate_ft(const double t,
		   const Vector<double>& v,
		   const double tolerance,
		   Vector<double>& result) const
  {
    result = 0;
  }

  void solve_ROW_stage_equation(const double t,
				const Vector<double>& v,
				const double alpha,
				const Vector<double>& y,
				const double tolerancs,
				Vector<double>& result) const
  {
//     // Jv=lambda*v -> (alpha*I-J)x=(alpha-lambda)x
//     result[0] = y[0] / (alpha - lambda_);
  }

  double period() const { return 17.0652165601579625588917206249; }

private:
  double lambda_;
};



int main()
{
  cout << "Testing one-step schemes for ODEs..." << endl;
  
  typedef Vector<double> V;

//   Dahlquist problem(5.0);
  Circle problem;

#if 0
  cout << "- checking consistency of the builtin one-step schemes:" << endl;

  cout << "* testing RK12:" << endl;
  ExplicitRungeKuttaScheme<V> rk12(ExplicitRungeKuttaScheme<V>::RK12);
  OneStepScheme<V>* scheme = &rk12;
//   Vector<double> temp(1), result(1), error_estimate(1);
  Vector<double> temp(2), result(2), error_estimate(2), exact(2);
  double err, olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }
  
//   cout << "* testing RK23:" << endl;
//   ExplicitRungeKuttaScheme<V> rk23(ExplicitRungeKuttaScheme<V>::RK23);
//   scheme = &rk23;
//   olderr = 0;
//   for (int expo = 0; expo <= 10; expo++) {
//     temp = problem.u0;
//     int N = 1<<expo;
//     double h = 1.0/N;
//     double err_est = 0;
//     for (int i = 1; i <= N; i++) {
//       scheme->increment(&problem, i*h, temp, h, result, error_estimate);
//       err_est = std::max(err_est, l2_norm(error_estimate));
//       temp = result;
//     }
//     problem.exact_solution(1.0, exact);
//     err = linfty_norm(result - exact);
//     if (expo > 0) {
//       cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
// 	   << ", max. of local error estimate " << err_est << endl;
//     }
//     olderr = err;
//   }
  
//   cout << "* testing Fehlberg34:" << endl;
//   ExplicitRungeKuttaScheme<V> fb34(ExplicitRungeKuttaScheme<V>::Fehlberg34);
//   scheme = &fb34;
//   olderr = 0;
//   for (int expo = 0; expo <= 10; expo++) {
//     temp = problem.u0;
//     int N = 1<<expo;
//     double h = 1.0/N;
//     double err_est = 0;
//     for (int i = 1; i <= N; i++) {
//       scheme->increment(&problem, i*h, temp, h, result, error_estimate);
//       err_est = std::max(err_est, l2_norm(error_estimate));
//       temp = result;
//     }
//     problem.exact_solution(1.0, exact);
//     err = linfty_norm(result - exact);
//     if (expo > 0) {
//       cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
// 	   << ", max. of local error estimate " << err_est << endl;
//     }
//     olderr = err;
//   }

//   cout << "* testing DoPri45:" << endl;
//   ExplicitRungeKuttaScheme<V> dopri45(ExplicitRungeKuttaScheme<V>::DoPri45);
//   scheme = &dopri45;
//   olderr = 0;
//   for (int expo = 0; expo <= 10; expo++) {
//     temp = problem.u0;
//     int N = 1<<expo;
//     double h = 1.0/N;
//     double err_est = 0;
//     for (int i = 1; i <= N; i++) {
//       scheme->increment(&problem, i*h, temp, h, result, error_estimate);
//       err_est = std::max(err_est, l2_norm(error_estimate));
//       temp = result;
//     }
//     problem.exact_solution(1.0, exact);
//     err = linfty_norm(result - exact);
//     if (expo > 0) {
//       cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
// 	   << ", max. of local error estimate " << err_est << endl;
//     }
//     olderr = err;
//   }

//   cout << "* testing DoPri78:" << endl;
//   ExplicitRungeKuttaScheme<V> dopri78(ExplicitRungeKuttaScheme<V>::DoPri78);
//   scheme = &dopri78;
//   olderr = 0;
//   for (int expo = 0; expo <= 10; expo++) {
//     temp = problem.u0;
//     int N = 1<<expo;
//     double h = 1.0/N;
//     double err_est = 0;
//     for (int i = 1; i <= N; i++) {
//       scheme->increment(&problem, i*h, temp, h, result, error_estimate);
//       err_est = std::max(err_est, l2_norm(error_estimate));
//       temp = result;
//     }
//     problem.exact_solution(1.0, exact);
//     err = linfty_norm(result - exact);
//     if (expo > 0) {
//       cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
// 	   << ", max. of local error estimate " << err_est << endl;
//     }
//     olderr = err;
//   }

  cout << "* testing ROS2:" << endl;
  ROWMethod<V> ros2(WMethod<V>::ROS2);
  scheme = &ros2;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing ROWDA3:" << endl;
  ROWMethod<V> rowda3(WMethod<V>::ROWDA3);
  scheme = &rowda3;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing ROS3:" << endl;
  ROWMethod<V> ros3(WMethod<V>::ROS3);
  scheme = &ros3;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing RODAS3:" << endl;
  ROWMethod<V> rodas3(WMethod<V>::RODAS3);
  scheme = &rodas3;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing GRK4T:" << endl;
  ROWMethod<V> grk4t(WMethod<V>::GRK4T);
  scheme = &grk4t;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing ROS3P:" << endl;
  ROWMethod<V> ros3p(WMethod<V>::ROS3P);
  scheme = &ros3p;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing ROS3Pw:" << endl;
  ROWMethod<V> ros3pw(WMethod<V>::ROS3Pw);
  scheme = &ros3pw;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing ROSI2P2:" << endl;
  ROWMethod<V> rosi2p2(WMethod<V>::ROSI2P2);
  scheme = &rosi2p2;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing RODASP:" << endl;
  ROWMethod<V> rodasp(WMethod<V>::RODASP);
  scheme = &rodasp;
  olderr = 0;
  for (int expo = 0; expo <= 10; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    problem.exact_solution(1.0, exact);
    err = linfty_norm(result - exact);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }
#endif

#if 0
  cout << "- checking adaptive solution of the Dahlquist test problem:" << endl;

  const double T = 1.0;
  const double q = 5.0;
  const double TOL = 1e-7;
  const double tau_max = 1.0;

  cout << "* TOL=" << TOL << endl;

  double errhelp;
  std::list<double>::const_iterator ti;
  IVPSolution<V> result_adaptive;

//   cout << "* testing RK12..." << endl;
//   ExplicitRungeKuttaScheme<V> rk12_adaptive(ExplicitRungeKuttaScheme<V>::RK12);
//   solve_IVP(&problem, &rk12_adaptive, T,
// 	    TOL, 0, q, tau_max, result_adaptive);

//   errhelp = 0;
//   ti = result_adaptive.t.begin();
//   for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
//        ui != result_adaptive.u.end(); ++ui, ++ti) {
//     cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
//   }
  
  cout << "* testing DoPri45..." << endl;
  ExplicitRungeKuttaScheme<V> dopri45_adaptive(ExplicitRungeKuttaScheme<V>::DoPri45);
  solve_IVP(&problem, &dopri45_adaptive, T,
	    TOL, 0, q, tau_max, result_adaptive);

  errhelp = 0;
  ti = result_adaptive.t.begin();
  for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
       ui != result_adaptive.u.end(); ++ui, ++ti) {
    cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
  }

  cout << "* testing DoPri78..." << endl;
  ExplicitRungeKuttaScheme<V> dopri78_adaptive(ExplicitRungeKuttaScheme<V>::DoPri78);
  solve_IVP(&problem, &dopri78_adaptive, T,
	    TOL, 0, q, tau_max, result_adaptive);

  errhelp = 0;
  ti = result_adaptive.t.begin();
  for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
       ui != result_adaptive.u.end(); ++ui, ++ti) {
    cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
  }

  
//   cout << "* testing ROS2..." << endl;
//   ROWMethod<V> ros2_adaptive(WMethod<V>::ROS2);
//   solve_IVP(&problem, &ros2_adaptive, T,
// 	    TOL, 0, q, tau_max, result_adaptive);
  
//   errhelp = 0;
//   ti = result_adaptive.t.begin();
//   for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
//        ui != result_adaptive.u.end(); ++ui, ++ti) {
//     cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
//   }
  
//   cout << "* testing ROS3..." << endl;
//   ROWMethod<V> ros3_adaptive(WMethod<V>::ROS3);
//   solve_IVP(&problem, &ros3_adaptive, T,
// 	    TOL, 0, q, tau_max, result_adaptive);
  
//   errhelp = 0;
//   ti = result_adaptive.t.begin();
//   for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
//        ui != result_adaptive.u.end(); ++ui, ++ti) {
//     cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
//   }
  
  cout << "* testing ROS3P..." << endl;
  ROWMethod<V> ros3p_adaptive(WMethod<V>::ROS3P);
  solve_IVP(&problem, &ros3p_adaptive, T,
	    TOL, 0, q, tau_max, result_adaptive);

  errhelp = 0;
  ti = result_adaptive.t.begin();
  for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
       ui != result_adaptive.u.end(); ++ui, ++ti) {
    cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
  }
#endif

#if 0
  cout << "- checking adaptive solution of the Arenstorf orbit problem:" << endl;

  Arenstorf problem2;

  const double T2 = problem2.period();
  const double q2 = 5.0;
  const double TOL2 = 1e-5;
  const double tau_max2 = 1.0;

  cout << "* TOL=" << TOL2 << endl;

  std::list<double>::const_iterator ti2;
  IVPSolution<V> result_adaptive2;

//   cout << "* testing RK12..." << endl;
//   ExplicitRungeKuttaScheme<V> rk12_adaptive2(ExplicitRungeKuttaScheme<V>::RK12);
//   solve_IVP(&problem2, &rk12_adaptive2, T2,
// 	    TOL2, 0, q2, tau_max2, result_adaptive2);
//   cout << "  ... done: " << result_adaptive2.u.size() << " time steps needed," << endl;
//   cout << "  initial value at T=0: "
//        << problem2.u0 << endl
//        << "  solution at  T=" << T2 << ": "
//        << *(result_adaptive2.u.rbegin())
//        << ", absolute error " << linfty_norm(problem2.u0 - *(result_adaptive2.u.rbegin())) << endl;
  
  cout << "* testing DoPri45..." << endl;
  ExplicitRungeKuttaScheme<V> dopri45_adaptive2(ExplicitRungeKuttaScheme<V>::DoPri45);
  solve_IVP(&problem2, &dopri45_adaptive2, T2,
 	    TOL2, TOL2, q2, tau_max2, result_adaptive2);
  cout << "  ... done: " << result_adaptive2.u.size() << " time steps needed," << endl;
  cout << "  initial value at T=0: "
       << problem2.u0 << endl
       << "  solution at  T=" << T2 << ": "
       << *(result_adaptive2.u.rbegin())
       << ", absolute error " << linfty_norm(problem2.u0 - *(result_adaptive2.u.rbegin())) << endl;
  
  cout << "* testing DoPri78..." << endl;
  ExplicitRungeKuttaScheme<V> dopri78_adaptive2(ExplicitRungeKuttaScheme<V>::DoPri78);
  solve_IVP(&problem2, &dopri78_adaptive2, T2,
	    TOL2, TOL2, q2, tau_max2, result_adaptive2);
  cout << "  ... done: " << result_adaptive2.u.size() << " time steps needed," << endl;
  cout << "  initial value at T=0: "
       << problem2.u0 << endl
       << "  solution at  T=" << T2 << ": "
       << *(result_adaptive2.u.rbegin())
       << ", absolute error " << linfty_norm(problem2.u0 - *(result_adaptive2.u.rbegin())) << endl;

//   cout << "* testing ROS2..." << endl;
//   ROWMethod<V> ros2_adaptive(WMethod<V>::ROS2);
//   solve_IVP(&problem, &ros2_adaptive, T,
// 	    TOL, 0, q, tau_max, result_adaptive);
  
//   errhelp = 0;
//   ti = result_adaptive.t.begin();
//   for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
//        ui != result_adaptive.u.end(); ++ui, ++ti) {
//     cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
//   }
  
//   cout << "* testing ROS3..." << endl;
//   ROWMethod<V> ros3_adaptive(WMethod<V>::ROS3);
//   solve_IVP(&problem, &ros3_adaptive, T,
// 	    TOL, 0, q, tau_max, result_adaptive);
  
//   errhelp = 0;
//   ti = result_adaptive.t.begin();
//   for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
//        ui != result_adaptive.u.end(); ++ui, ++ti) {
//     cout << "  absolute error at t=" << *ti << ": " << fabs((*ui)[0] - problem.exact_solution(*ti)) << endl;
//   }

//   cout << "* testing ROS3P..." << endl;
//   ROWMethod<V> ros3p_adaptive2(WMethod<V>::ROS3P);
//   solve_IVP(&problem2, &ros3p_adaptive2, T2,
// 	    TOL2, TOL2, q2, tau_max2, result_adaptive2);
//   cout << "  ... done: " << result_adaptive2.u.size() << " time steps needed," << endl;
//   cout << "  initial value at T=0: "
//        << problem2.u0 << endl
//        << "  solution at  T=" << T2 << ": "
//        << *(result_adaptive2.u.rbegin())
//        << ", absolute error " << linfty_norm(problem2.u0 - *(result_adaptive2.u.rbegin())) << endl;
  
#endif

#if 1
  cout << "- checking the transform_coefficients() routine:" << endl;
  LowerTriangularMatrix<double> Alpha, Gamma, A, C;
  Vector<double> b, bhat, m, e, alpha_vector, gamma_vector;
  unsigned int s;

  cout << "* ROS2:" << endl;
  s = 2;
  Alpha.resize(2, 2);
  cout << "Alpha=" << endl << Alpha;

  ROWMethod<V>::transform_coefficients(Alpha, Gamma, b, bhat, A, C, m, e, alpha_vector, gamma_vector);
#endif

  return 0;
}

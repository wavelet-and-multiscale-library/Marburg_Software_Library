#include <cmath>
#include <iostream>
#include <algebra/vector.h>
#include <numerics/ivp.h>
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
    // Jv=lambda*v -> (I-alpha*J)x=(1-alpha*lambda)x
    result[0] = y[0] / (1 - alpha*lambda_);
  }

private:
  double lambda_;
};


int main()
{
  cout << "Testing one-step schemes for ODEs..." << endl;
  
  typedef Vector<double> V;

  Dahlquist problem;

  cout << "* testing RK12:" << endl;
  ExplicitRungeKuttaScheme<V> rk12(ExplicitRungeKuttaScheme<V>::RK12);
  OneStepScheme<V>* scheme = &rk12;
  Vector<double> temp(1), result(1), error_estimate(1);
  double err, olderr = 0;
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }
  
  cout << "* testing RK23:" << endl;
  ExplicitRungeKuttaScheme<V> rk23(ExplicitRungeKuttaScheme<V>::RK23);
  scheme = &rk23;
  olderr = 0;
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }
  
  cout << "* testing Fehlberg34:" << endl;
  ExplicitRungeKuttaScheme<V> fb34(ExplicitRungeKuttaScheme<V>::Fehlberg34);
  scheme = &fb34;
  olderr = 0;
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing DoPri45:" << endl;
  ExplicitRungeKuttaScheme<V> dopri45(ExplicitRungeKuttaScheme<V>::DoPri45);
  scheme = &dopri45;
  olderr = 0;
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing DoPri78:" << endl;
  ExplicitRungeKuttaScheme<V> dopri78(ExplicitRungeKuttaScheme<V>::DoPri78);
  scheme = &dopri78;
  olderr = 0;
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  cout << "* testing ROS2:" << endl;
  ROWMethod<V> ros2(WMethod<V>::ROS2);
  scheme = &ros2;
  olderr = 0;
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    err = fabs(result[0] - M_E);
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
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    double err_est = 0;
    for (int i = 1; i <= N; i++) {
      scheme->increment(&problem, i*h, temp, h, result, error_estimate);
      err_est = std::max(err_est, l2_norm(error_estimate));
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2
	   << ", max. of local error estimate " << err_est << endl;
    }
    olderr = err;
  }

  return 0;
}

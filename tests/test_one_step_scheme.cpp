#include <iostream>
#include <numerics/ivp.h>
#include <numerics/one_step_scheme.h>
#include <numerics/runge_kutta.h>
#include <algebra/vector.h>

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

  void apply_f(const double t, const Vector<double>& v, const double tolerance,
	       Vector<double>& result) const
  {
    result[0] = lambda_ * v[0];
  }

private:
  double lambda_;
};

//   cout << "* testing ROS2:" << endl;
//   Rosenbrock<Point<1>, IVP<1> > R2(Rosenbrock<Point<1>, IVP<1> >::ROS2);
//   for (int expo = 0; expo <= 6; expo++)
//     {
//       temp = problem.u0;
//       int N = 1<<expo;
//       double h = 1.0/N;
//       for (int i = 1; i <= N; i++)
// 	{
// 	  R2.increment(problem, i*h, temp, h, result);
// 	  temp = result;
// 	}
//       err = fabs(result[0] - M_E);
//       if (expo > 0)
// 	{
// 	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
// 	}
//       olderr = err;
//     }

//   cout << "* testing ROS3:" << endl;
//   Rosenbrock<Point<1>, IVP<1> > R3(Rosenbrock<Point<1>, IVP<1> >::ROS3);
//   for (int expo = 0; expo <= 6; expo++)
//     {
//       temp = problem.u0;
//       int N = 1<<expo;
//       double h = 1.0/N;
//       for (int i = 1; i <= N; i++)
// 	{
// 	  R3.increment(problem, i*h, temp, h, result);
// 	  temp = result;
// 	}
//       err = fabs(result[0] - M_E);
//       if (expo > 0)
// 	{
// 	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
// 	}
//       olderr = err;
//     }

//   cout << "* testing ROWDA3:" << endl;
//   Rosenbrock<Point<1>, IVP<1> > R4(Rosenbrock<Point<1>, IVP<1> >::ROWDA3);
//   for (int expo = 0; expo <= 6; expo++)
//     {
//       temp = problem.u0;
//       int N = 1<<expo;
//       double h = 1.0/N;
//       for (int i = 1; i <= N; i++)
// 	{
// 	  R4.increment(problem, i*h, temp, h, result);
// 	  temp = result;
// 	}
//       err = fabs(result[0] - M_E);
//       if (expo > 0)
// 	{
// 	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
// 	}
//       olderr = err;
//     }

//   cout << "* testing RODAS3:" << endl;
//   Rosenbrock<Point<1>, IVP<1> > R5(Rosenbrock<Point<1>, IVP<1> >::RODAS3);
//   for (int expo = 0; expo <= 6; expo++)
//     {
//       temp = problem.u0;
//       int N = 1<<expo;
//       double h = 1.0/N;
//       for (int i = 1; i <= N; i++)
// 	{
// 	  R5.increment(problem, i*h, temp, h, result);
// 	  temp = result;
// 	}
//       err = fabs(result[0] - M_E);
//       if (expo > 0)
// 	{
// 	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
// 	}
//       olderr = err;
//     }

int main()
{
  cout << "Testing one-step schemes for ODEs..." << endl;
  
  typedef Vector<double> V;

  Dahlquist problem;

  cout << "* testing RK12:" << endl;
  ExplicitRungeKuttaScheme<V> rk12(ExplicitRungeKuttaScheme<V>::RK12);
  ExplicitRungeKuttaScheme<V>* scheme = &rk12;
  Vector<double> temp(1), result(1), error_estimate(1);
  double err, olderr = 0;
  for (int expo = 0; expo <= 6; expo++) {
    temp = problem.u0;
    int N = 1<<expo;
    double h = 1.0/N;
    for (int i = 1; i <= N; i++) {
      scheme->increment(problem, i*h, temp, h, result, error_estimate);
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
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
    for (int i = 1; i <= N; i++) {
      scheme->increment(problem, i*h, temp, h, result, error_estimate);
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
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
    for (int i = 1; i <= N; i++) {
      scheme->increment(problem, i*h, temp, h, result, error_estimate);
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
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
    for (int i = 1; i <= N; i++) {
      scheme->increment(problem, i*h, temp, h, result, error_estimate);
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
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
    for (int i = 1; i <= N; i++) {
      scheme->increment(problem, i*h, temp, h, result, error_estimate);
      temp = result;
    }
    err = fabs(result[0] - M_E);
    if (expo > 0) {
      cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
    }
    olderr = err;
  }

  return 0;
}

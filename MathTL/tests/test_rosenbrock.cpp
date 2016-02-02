#include <iostream>
#include <geometry/point.h>
#include <numerics/ivp.h>
#include <numerics/rosenbrock.h>

using std::cout;
using std::endl;
using namespace MathTL;

/*
  The Dahlquist standard problem
    u' = lambda*u, u(0) = 1
  
  So f(t,u)=lambda*u, f_t(t,u)=0 and f_u(t,u)=lambda*I.
 */
class Dahlquist
  : public IVP<1>
{
public:
  Dahlquist(const double lambda = 1)
    : lambda_(lambda)
  {
    u0[0] = 1;
  }
  
  void apply_f(const double t, const Point<1>& v,
	       Point<1>& result) const
  {
    result[0] = lambda_ * v[0];
  }

  void apply_ft(const double t, const Point<1>& v,
		Point<1>& result) const
  {
    result[0] = 0; // no t-dependence
  }

  void solve_jacobian(const double t, const Point<1>& v, const double tau,
		      Point<1>& result) const
  {
    result[0] = v[0] / (1 - lambda_*tau);
  }

private:
  double lambda_;
};

int main()
{
  cout << "Testing MathTL::Rosenbrock ..." << endl;

  Dahlquist problem;

  cout << "* testing the linear-implicit Euler:" << endl;
  Rosenbrock<Point<1>, IVP<1> > R1;
  Point<1> temp, result;
  double err, olderr = 0;
  for (int expo = 0; expo <= 6; expo++)
    {
      temp = problem.u0;
      int N = 1<<expo;
      double h = 1.0/N;
      for (int i = 1; i <= N; i++)
	{
	  R1.increment(problem, i*h, temp, h, result);
	  temp = result;
	}
      err = fabs(result[0] - M_E);
      if (expo > 0)
	{
	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
	}
      olderr = err;
    }

  cout << "* testing ROS2:" << endl;
  Rosenbrock<Point<1>, IVP<1> > R2(Rosenbrock<Point<1>, IVP<1> >::ROS2);
  for (int expo = 0; expo <= 6; expo++)
    {
      temp = problem.u0;
      int N = 1<<expo;
      double h = 1.0/N;
      for (int i = 1; i <= N; i++)
	{
	  R2.increment(problem, i*h, temp, h, result);
	  temp = result;
	}
      err = fabs(result[0] - M_E);
      if (expo > 0)
	{
	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
	}
      olderr = err;
    }

  cout << "* testing ROS3:" << endl;
  Rosenbrock<Point<1>, IVP<1> > R3(Rosenbrock<Point<1>, IVP<1> >::ROS3);
  for (int expo = 0; expo <= 6; expo++)
    {
      temp = problem.u0;
      int N = 1<<expo;
      double h = 1.0/N;
      for (int i = 1; i <= N; i++)
	{
	  R3.increment(problem, i*h, temp, h, result);
	  temp = result;
	}
      err = fabs(result[0] - M_E);
      if (expo > 0)
	{
	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
	}
      olderr = err;
    }

  cout << "* testing ROWDA3:" << endl;
  Rosenbrock<Point<1>, IVP<1> > R4(Rosenbrock<Point<1>, IVP<1> >::ROWDA3);
  for (int expo = 0; expo <= 6; expo++)
    {
      temp = problem.u0;
      int N = 1<<expo;
      double h = 1.0/N;
      for (int i = 1; i <= N; i++)
	{
	  R4.increment(problem, i*h, temp, h, result);
	  temp = result;
	}
      err = fabs(result[0] - M_E);
      if (expo > 0)
	{
	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
	}
      olderr = err;
    }

  cout << "* testing RODAS3:" << endl;
  Rosenbrock<Point<1>, IVP<1> > R5(Rosenbrock<Point<1>, IVP<1> >::RODAS3);
  for (int expo = 0; expo <= 6; expo++)
    {
      temp = problem.u0;
      int N = 1<<expo;
      double h = 1.0/N;
      for (int i = 1; i <= N; i++)
	{
	  R5.increment(problem, i*h, temp, h, result);
	  temp = result;
	}
      err = fabs(result[0] - M_E);
      if (expo > 0)
	{
	  cout << "h=" << h << ", error " << err << ", p approx. " << log(olderr/err)/M_LN2 << endl;
	}
      olderr = err;
    }

  return 0;
}

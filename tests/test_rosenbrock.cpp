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
  Rosenbrock<Point<1>, IVP<1> > R1;

  Point<1> result;
  cout << problem.u0 << endl;
  R1.increment(problem, 0, problem.u0, 0.1, result);
  cout << result << endl;

  return 0;
}

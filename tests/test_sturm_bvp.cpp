#include <iostream>
#include <numerics/bvp.h>
#include <numerics/sturm_bvp.h>

using std::cout;
using std::endl;
using namespace MathTL;

/*
  -y''=1, y(0)=y(1)=0
 */
class TestProblem
  : public SturmBVP
{
public:
  double p(const double t) const { return 1; }
  double p_prime(const double t) const { return 0; }
  double q(const double t) const { return 0; }
  double g(const double t) const { return 1; }
  double alpha0() const { return 1; }
  double alpha1() const { return 0; }
  double beta0() const { return 1; }
  double beta1() const { return 0; }
};

int main()
{
  cout << "Testing MathTL::SturmBVP ..." << endl;

  TestProblem T;

  return 0;
}

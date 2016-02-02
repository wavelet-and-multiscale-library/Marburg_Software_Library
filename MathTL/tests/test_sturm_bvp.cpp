#include <iostream>
#include <numerics/bvp.h>
#include <numerics/sturm_bvp.h>

using std::cout;
using std::endl;
using namespace MathTL;

/*
  a test problem
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

/*
  the same problem
    -y''=1, y(0)=y(1)=0
  now as SimpleSturmBVP
 */
class TestProblem2
  : public SimpleSturmBVP
{
public:
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
  double p(const double t) const { return 1; }
  double p_prime(const double t) const { return 0; }
  double q(const double t) const { return 0; }
  double g(const double t) const { return 1; }
};

int main()
{
  cout << "Testing MathTL::SturmBVP ..." << endl;

  TestProblem T;
  TestProblem2 T2;

  return 0;
}

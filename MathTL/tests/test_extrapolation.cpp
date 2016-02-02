#include <iostream>
#include <numerics/extrapolation.h>
#include <numerics/quadrature.h>
#include <utils/function.h>

using std::cout;
using std::endl;
using namespace MathTL;

// example 9.26 in Deuflhard/Hohmann, Numerische Mathematik I
// quadrature of int_{-1}^1 1/(10^{-4}+x^2) dx

class TestFunction : public Function<1>
{
public:
  TestFunction() : Function<1>(1) {}
  virtual ~TestFunction(){}
  double value(const Point<1>& p,
	       const unsigned int component = 0) const
  {
    return 1.0/((1e-4)+p[0]*p[0]);
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

class Romberg
{
public:
  Romberg(const Function<1>& f) : f_(f) {}
  void approximate(const unsigned int n, double& r) const
  {
    CompositeRule<1> composite(trapez, n);
    r = composite.integrate(f_, Point<1>(-1.0), Point<1>(1.0));
  }
  
protected:
  const Function<1>& f_;
  TrapezoidalRule trapez;
};

int main()
{
  cout << "Testing extrapolation methods ..." << endl;

  RombergSequence r;
  cout << "- the first 10 entries of the Romberg sequence:" << endl;
    for (unsigned int k = 1; k <= 10; k++)
    cout << r.n(k) << " ";
  cout << endl;

  BulirschSequence b;
  cout << "- the first 10 entries of the Bulirsch sequence:" << endl;
    for (unsigned int k = 1; k <= 10; k++)
    cout << b.n(k) << " ";
  cout << endl;

  TestFunction f;
  Romberg T(f);
  ExtrapolationTable<Romberg> E(T, 13, 2);

  cout << "- the extrapolation table for the Deuflhard/Hohmann example:" << endl
       << E.table();

  return 0;
}

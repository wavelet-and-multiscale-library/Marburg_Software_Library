#include <iostream>
#include <algebra/polynomial.h>

using namespace std;
using namespace MathTL;

int main(int, char **)
{
  Polynomial<double> p, q;

  cout << "Testing class Polynomial..." << endl;

  cout << "The zero polynomial:" << endl << p << endl << "... has degree " << p.degree() << endl;
  
  p.set_coefficient(0, 23);
  cout << "A constant polynomial:" << endl << p << endl;

  p.set_coefficient(2, -2);
  cout << "A non-constant polynomial p:" << endl << p << endl;

  double x = 1.2;
  cout << "Evaluating p at x=" << x << ":" << endl << p.value(x)
       << " (with full Horner scheme: " << p.value(x, 0) << ")"
       << endl;

  q.set_coefficient(0, -20);
  q.set_coefficient(1, 4);
  cout << "Another polynomial q:" << endl << q << endl;

  p.scale(2.0);
  cout << "Substitution of x=2y in p:" << endl << p << endl;
  q.scale(2.0);
  cout << "Substitution of x=2y in q:" << endl << q << endl;

  p.shift(1.0);
  cout << "Shift x=y+1 in p:" << endl << p << endl;

  q.chain(p);
  cout << "q circ p:" << endl << q << endl;

  q *= -2;
  cout << "q*=-2:" << endl << q << endl;

  p += q;
  cout << "p+=q:" << endl << p << endl;

  p -= q;
  cout << "p-=q:" << endl << p << endl;

  p *= q;
  cout << "p*=q:" << endl << p << endl;

  cout << "p+q and p:" << endl << p+q << endl << p << endl;

  cout << "p-q and p:" << endl << p-q << endl << p << endl;

  cout << "p*q and p:" << endl << p*q << endl << p << endl;
  
  cout << "-p and p:" << endl << -p << endl << p << endl;

  cout << "-2*p and p:" << endl << (-2.0)*p << endl << p << endl;

  cout << "Some powers of p:" << endl;
  for (unsigned int k(0); k <= 2; k++)
    cout << "p^" << k << ":" << endl << p.power(k) << endl;

  cout << "derivative of p:" << endl << p.differentiate() << endl;
  cout << "antiderivative of p: " << endl << p.integrate() << endl;

  Polynomial<double> r;
  r.set_coefficient(2, 1);
  double a = 0.5;
  double b = 1.4;
  cout << "Integrating the polynomial" << endl << r << endl << "over the interval ["
       << a << "," << b << "]:" << endl << r.integrate(a, b) << endl;

  cout << "And now with a quadrature rule:" << endl
       << r.integrate(a, b, true) << endl;

  cout << "Polynomial p again:" << endl << p << endl;
  cout << "p'(3.1415)=" << p.value(3.1415, 1) << "="
       << p.differentiate().value(3.1415) << endl;
  cout << "p''(3.1415)=" << p.value(3.1415, 2) << "="
       << p.differentiate().differentiate().value(3.1415) << endl;
}

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
  cout << "Evaluating p at x=" << x << ":" << endl << p(x) << endl;

  q.set_coefficient(0, -20);
  q.set_coefficient(1, 4);
  cout << "Another polynomial q:" << endl << q << endl;

//   cout << "p circ q:" << endl << p.substituteInto(q) << endl;

//   q *= -2;
//   cout << "q*=-2:" << endl << q << endl;

//   p += q;
//   cout << "p+=q:" << endl << p << endl;

//   p -= q;
//   cout << "p-=q:" << endl << p << endl;

//   p *= q;
//   cout << "p*=q:" << endl << p << endl;

//   cout << "p+q:" << endl << p+q << endl << p << endl;

//   cout << "p-q:" << endl << p-q << endl << p << endl;

//   cout << "p*q:" << endl << p*q << endl << p << endl;
  
//   cout << "-p:" << endl << -p << endl << p << endl;

//   cout << "-2*p:" << endl << (-2.0)*p << endl << p << endl;

//   cout << "Some powers of p:" << endl;
//   for (int k(0); k <= 2; k++)
//     cout << "p^" << k << ":" << endl << p.power(k) << endl;

//   cout << "derivative of p:" << endl << p.sdifferentiate() << endl;
//   cout << "antiderivative of p: " << endl << p.sintegrate() << endl;

//   Polynomial r;
//   r.setCoeff(2, 1);
//   double a = 0.5;
//   double b = 1.4;
//   cout << endl << "Integrating the polynomial" << endl << r << endl << "over the interval ["
//        << a << "," << b << "]:" << endl << r.integrate(a, b) << endl;

//   cout << endl << "And now with quadrature rule:" << endl
//        << r.integrate(a, b, true) << endl;

//   cout << endl;

//   cout << "Polynomial p again:" << endl << p << endl;
//   cout << "p'(3.14)=" << p.derivative(3.14) << "=" << (p.sdifferentiate())(3.14) << endl;
}

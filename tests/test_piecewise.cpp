#include <iostream>
#include <cmath>
#include <algebra/piecewise.h>
#include <algebra/polynomial.h>
#include <utils/array1d.h>
#include <numerics/gram_schmidt.h>

using namespace std;
using namespace MathTL;

int main()
{
  cout << "Testing Piecewise ..." << endl;

  Piecewise<double> p, r;
  Polynomial<double> q;
  q = 2;
  p.set_local_expansion(-1, q);

  q = 4;
  p.set_local_expansion(2, -q);
  
  r.set_local_expansion(1, (-2.0)*q);
  r.set_local_expansion(2, q);

  cout << "A simple piecewise function p:" << endl << p << endl;
  cout << "Another one r:" << endl << r << endl;

  int j = 1;
  p.split_me(j);
  cout << "Splitting p at granularity j=" << j << ":" << endl << p << endl;

  int k = 1;
  cout << "r shifted by k=" << k << ":" << endl << r.shift(k) << "r remained unchanged:" << endl << r << endl;

  int k1=3;
  int k2=5;
  cout << "p clipped to [" << ldexp(1.0, -p.get_granularity())*k1
       << "," << ldexp(1.0, -p.get_granularity())*k2 << "]:" << endl
       << p.clip(k1, k2) << "p remained unchanged:" << endl << p << endl;

  cout << "p+r:" << endl << p+r << "p remained unchanged:" << endl << p << endl;
  cout << "p*r:" << endl << p*r << "p remained unchanged:" << endl << p << endl;

  Piecewise<double> s(2);
  q.set_coefficient(0, -1);
  q.set_coefficient(2, 1);
  s.set_local_expansion(0, q);
  q = 1;
  q.set_coefficient(1, -2);
  s.set_local_expansion(2, q);
  cout << "A more subtle piecewise function s:" << endl << s << endl;
  cout << "s':" << endl << s.differentiate() << endl;
  cout << "int(s(x),x=-infinity..infinity):" << endl << s.integrate() << endl;

  cout << "s'(0.14) symbolically/direct:" << endl << s.differentiate()(0.14)
       << "/" << s.derivative(0.14) << endl;


//   double x=2.7;
//   cout << "At x=" << x << ", it takes the value" << endl << p(x) << endl;

//   p *= 2.0;
//   cout << "p*=2:" << endl << p << endl;

//   p *= q;
//   cout << "p*=q:" << endl << p << endl;

//   p += q;
//   cout << "p+=q:" << endl << p << endl;

//   p *= p;
//   cout << "p*=p:" << endl << p << endl;

//   cout << "p+r:" << endl << p+r << endl << p << endl;
//   cout << "2.5*p-r:" << endl << 2.5*p-r << endl << p << endl;
//   cout << "p*r:" << endl << p*r << endl << p << endl;

//   int k=2;
//   cout << "p shifted by k=" << k << ":" << endl << p.shift(k) << endl << p << endl;

  // monomial basis
  Polynomial<double> x0;
  x0.set_coefficient(0,1);
  Polynomial<double> x1;
  x1.set_coefficient(1,1);
  Polynomial<double> x2;
  x2.set_coefficient(2,1);

  // ... on [-1,1]
  Piecewise<double> y0, y1, y2;
  y0.set_local_expansion(-1, x0);
  y0.set_local_expansion(0, x0);
  y1.set_local_expansion(-1, x1);
  y1.set_local_expansion(0, x1);
  y2.set_local_expansion(-1, x2);
  y2.set_local_expansion(0, x2);
  cout << "monomial basis on [-1,1] ..." << endl;
  cout << "y0 = " << y0;
  cout << "y1 = " << y1;
  cout << "y2 = " << y2;

  Array1D< Piecewise<double> > b(3);
  b[0] = y0;
  b[1] = y1;
  b[2] = y2;
  gramSchmidtProcess(b);
  cout << "Gram-Schmidt process ..." << endl;
  cout << b << endl;

  return 0;
}

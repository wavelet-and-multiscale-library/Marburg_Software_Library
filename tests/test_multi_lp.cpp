#include <iostream>
#include <MathTL.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MultivariateLaurentPolynomial<>..." << endl;

  cout << "- the zero multivariate Laurent polynomial:" << endl;
  typedef MultivariateLaurentPolynomial<double, 2> LPoly;
  LPoly p;
  cout << p << endl;
  
  cout << "- a constant multivariate Laurent polynomial:" << endl;
  LPoly r(42);
  cout << r << endl;

  cout << "- checking assignment operator (q=r):" << endl;
  LPoly q = r;
  cout << q << endl;
  
  cout << "- adding a nontrivial power:" << endl;
  MultiIndex<int, 2> index;
  index[0] = 1; index[1] = 2;
  p.set_coefficient(index, 23);
  cout << p << endl;
  
  cout << "- adding another power:" << endl;
  index[0] = -3; index[1] = 3;
  p.set_coefficient(index, 4);
  cout << p << endl;
  
  cout << "- evaluating p at some (nontrivial) points:" << endl;
  Point<2> x(1, 1);
  cout << "p(" << x << ")=" << p.value(x) << endl;
  x = Point<2>(2, -3);
  cout << "p(" << x << ")=" << p.value(x) << endl;

  cout << "- algebraic functionality:" << endl;
  q += p;
//   q.dump();
  cout << "q+=p: " << q << endl;

  q -= p;
  q -= p;
  cout << "q-=p; q-=p: " << q << endl;
  q = p;
  q *= 3;
  cout << "q=p; q*=3: " << q << endl;
  cout << "q+p: " << q+p << endl;
  cout << "q-p: " << q-p << endl;
  cout << "-q: " << -q << endl;
  cout << "3*q: " << 3. * q << endl;

//   cout << "p*q: " << p*q << endl;
//   cout << "(1+z)^2: ";
//   LaurentPolynomial<double> s(1);
//   s.set_coefficient(1, 1);
//   cout << s.power(2) << endl;
//   cout << "(1+z)^3: " << s.power(3) << endl;

//   LaurentPolynomial<double> t;
//   t.set_coefficient(-1, 43);
//   t.set_coefficient(0, 3);
//   t.set_coefficient(1, 3);
//   t.set_coefficient(2, 1);
//   cout << "a new Laurent polynomial t:" << endl << t << endl;
//   LaurentPolynomial<double> divisor, dividend, remainder;
//   divisor.set_coefficient(0, 1);
//   divisor.set_coefficient(1, 1);
//   t.divide(divisor, dividend, remainder);
//   cout << "division of t by " << divisor << " yields dividend" << endl
//        << dividend << endl
//        << "and remainder" << endl
//        << remainder << endl;
  
  return 0;
}

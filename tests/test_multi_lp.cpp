#include <iostream>
#include <MathTL.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing MultiLaurentPolynomial<>..." << endl;

  cout << "- the zero multivariate Laurent polynomial:" << endl;
  typedef MultiLaurentPolynomial<double, 2> LPoly;
  LPoly p;
//   cout << p << endl;
//   cout << "  (it has degree " << p.degree() << ")" << endl;
  
//   cout << "- a constant polynomial r\\in\\mathbb L[x], x\\in\\mathbb R:" << endl;
//   LaurentPolynomial<double> r(42);
//   cout << r << endl;
//   cout << "  (it has degree " << r.degree() << ")" << endl;

//   cout << "- checking assignment operator (q=r):" << endl;
//   LaurentPolynomial<double> q = r;
//   cout << q << endl;
  
//   cout << "- adding a nontrivial power:" << endl;
//   p.set_coefficient(3, 23);
//   cout << p << endl;
  
//   cout << "- adding a negative power:" << endl;
//   p.set_coefficient(-2, 4);
//   cout << p << endl;
//   cout << "  (p now has degree " << p.degree() << ")" << endl;
  
//   cout << "- evaluating p at some (nontrivial) points:" << endl;
//   double x(1);
//   cout << "p(" << x << ")=" << p.value(x) << endl;
//   x = 2;
//   cout << "p(" << x << ")=" << p.value(x) << endl;
//   x = -0.3;
//   cout << "p(" << x << ")=" << p.value(x) << endl;

//   cout << "- algebraic functionality:" << endl;
//   q += p;
//   cout << "q+=p: " << q << endl;
//   q -= p;
//   q -= p;
//   cout << "q-=p; q-=p: " << q << endl;
//   q = p;
//   q *= 3;
//   cout << "q=p; q*=3: " << q << endl;
//   cout << "q+p: " << q+p << endl;
//   cout << "q-p: " << q-p << endl;
//   cout << "-q: " << -q << endl;
//   cout << "3*q: " << 3.*q << endl;
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

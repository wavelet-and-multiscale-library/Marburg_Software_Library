#include <iostream>
#include <MathTL.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing LaurentPolynomial<>..." << endl;

  cout << "- the zero polynomial p\\in\\mathbb L[x], x\\in\\mathbb R:" << endl;
  LaurentPolynomial<double> p;
  cout << p << endl;
  
  cout << "- a constant polynomial r\\in\\mathbb L[x], x\\in\\mathbb R:" << endl;
  LaurentPolynomial<double> r(42);
  cout << r << endl;

  cout << "- checking assignment operator (q=r):" << endl;
  LaurentPolynomial<double> q = r;
  cout << q << endl;
  
  cout << "- adding a nontrivial power:" << endl;
  p.set_coefficient(3, 23);
  cout << p << endl;
  
  cout << "- adding a negative power:" << endl;
  p.set_coefficient(-2, 4);
  cout << p << endl;
  
//   cout << "- evaluating p at some (nontrivial) points:" << endl;
//   double x(1);
//   cout << "p(" << x << ")=" << p.value(x) << endl;
//   x = 2;
//   cout << "p(" << x << ")=" << p.value(x) << endl;
//   x = -0.3;
//   cout << "p(" << x << ")=" << p.value(x) << endl;

  cout << "- algebraic functionality:" << endl;
  q += p;
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
  cout << "3*q: " << 3.*q << endl;
//   cout << "p*q: " << p*q << endl;
//   cout << "(1+z)^2: ";
//   LaurentPolynomial<double> r;
//   r.setCoefficient(0, 1);
//   r.setCoefficient(1, 1);
//   cout << power(r, 2) << endl;

  return 0;
}

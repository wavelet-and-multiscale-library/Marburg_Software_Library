#include <iostream>
#include <cmath>
#include <numerics/ortho_poly.h>
#include <algebra/vector.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing OrthogonalPolynomial class..." << endl;

  double x;

  x = 0.4321;
  cout << "- evaluating monomials at x=" << x << ":" << endl;
  Monomial mon;
  for (int n(0); n <= 10; n++)
    cout << "  * T_" << n << "(" << x << ")=" << mon(n, x)
 	 << " (error: " << fabs(mon(n,x)-pow(x,n)) << ")" << endl;
  
  x = 0.1234;
  cout << "- evaluating Chebyshev polynomials at x=" << x << ":" << endl;
  ChebyshevPolynomial cheb;
  for (int n(0); n <= 6; n++)
    {
      double exact(0), approx(cheb(n,x));
      switch (n)
	{
	case 0: exact = 1.0; break;
	case 1: exact = x; break;
	case 2: exact = x*x - 0.5; break;
	case 3: exact = x*x*x - 3.0/4.0*x; break;
	case 4: exact = x*x*x*x - x*x + 1.0/8.0; break;
	case 5: exact = x*x*x*x*x - 20.0/16.0*x*x*x + 5.0/16.0*x; break;
	case 6: exact = x*x*x*x*x*x - 48.0/32.0*x*x*x*x + 18.0/32.0*x*x - 1.0/32.0; break;
	}
      cout << "  * T_" << n << "(" << x << ")=" << approx
	   << " (error: " << fabs(approx-exact) << ")" << endl;
    }

  x = 0.99;
  cout << "- evaluating Legendre polynomials at x=" << x << ":" << endl;
  LegendrePolynomial leg;
  for (int n(0); n <= 4; n++)
    {
      double exact(0), approx(leg(n,x));
      switch (n)
	{
	case 0: exact = 1.0; break;
	case 1: exact = x; break;
	case 2: exact = x*x - 1.0/3.0; break;
	case 3: exact = x*x*x - 3.0/5.0*x; break;
	case 4: exact = x*x*x*x - 6.0/7.0*x*x + 3.0/35.0; break;
	}
      cout << "  * T_" << n << "(" << x << ")=" << approx
	   << " (error: " << fabs(approx-exact) << ")" << endl;
    }

  Vector<double> coeffs(3), coeffs2(1);
  coeffs(0) = -18.0;
  coeffs(1) = 4.0;
  coeffs(2) = 1.0; // x^2+4*x-18
  x = 6.0;
  cout << "- forward summation of a polynomial: "
       << Monomial().forwardSummation(coeffs, x) << endl;
  cout << "- adjoint summation of a polynomial (Horner scheme): "
       << Monomial().adjointSummation(coeffs, x) << endl;
  coeffs2(0) = 23.0;
  cout << "- forward summation of another polynomial: "
       << Monomial().forwardSummation(coeffs2, x) << endl;
  cout << "- adjoint summation of another (Horner scheme): "
       << Monomial().adjointSummation(coeffs2, x) << endl;
  
  cout << "- test assembly of normalized orthogonal polynomials:" << endl;
  cout << "  * some monomials:" << endl;
  for (int n(0); n <= 4; n++)
    cout << Monomial().assemble(n) << endl;
  cout << "  * some Chebyshev polynomials:" << endl;
  for (int n(0); n <= 4; n++)
    cout << ChebyshevPolynomial().assemble(n) << endl;
  cout << "  * some Legendre polynomials:" << endl;
  for (int n(0); n <= 4; n++)
    cout << LegendrePolynomial().assemble(n) << endl;

  return 0;
}

#include <iostream>
#include <cmath>
#include <numerics/recursion.h>
#include <algebra/vector.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing ThreeTermRecursion class..." << endl;
  double x;

  x = 0.4321;
  cout << "- monomial recursion for x=" << x << ":" << endl;
  MonomialRecursion monrec(x);
  for (int n(0); n <= 10; n++)
    cout << "  * T_" << n << "(" << x << ")=" << monrec(n)
 	 << " (error: " << fabs(monrec(n)-pow(x,n)) << ")" << endl;

  x = 0.1;
  cout << "- Chebyshev recursion for x=" << x << ":" << endl;
  ChebyshevRecursion cr(x);
  for (int n(0); n <= 6; n++)
    {
      double exact(0), approx(cr(n));
      switch (n)
	{
	case 0: exact = 1.0; break;
	case 1: exact = x; break;
	case 2: exact = 2.0*x*x - 1.0; break;
	case 3: exact = 4.0*x*x*x - 3.0*x; break;
	case 4: exact = 8.0*x*x*x*x - 8.0*x*x + 1.0; break;
	case 5: exact = 16.0*x*x*x*x*x - 20.0*x*x*x + 5.0*x; break;
	case 6: exact = 32.0*x*x*x*x*x*x - 48.0*x*x*x*x + 18.0*x*x - 1.0; break;
	}
      cout << "  * T_" << n << "(" << x << ")=" << approx
	   << " (error: " << fabs(approx-exact) << ")" << endl;
    }

  x = 0.923;
  cout << "- checking Deuflhard/Hohmann exercise 6.1 for T_31: "
       << ChebyshevRecursion(x)(31)
       << " (error: " << fabs(ChebyshevRecursion(x)(31)-0.948715916161)
       << ")" << endl;

  x = 0.1234;
  cout << "- cosine recursion for x=" << x << ":" << endl;
  CosineRecursion cosrec(x);
  for (int n(0); n <= 10; n++)
    cout << "  * T_" << n << "(" << x << ")=" << cosrec(n)
	 << " (error: " << fabs(cosrec(n)-cos(n*x)) << ")" << endl;

  x = 0.1234;
  cout << "- sine recursion for x=" << x << ":" << endl;
  SineRecursion sinrec(x);
  for (int n(0); n <= 10; n++)
    cout << "  * T_" << n << "(" << x << ")=" << sinrec(n)
 	 << " (error: " << fabs(sinrec(n)-sin(n*x)) << ")" << endl;

  x = 0.1;
  const double a(1.0), b(1.0);
  cout << "- Jacobi recursion for x=" << x << ":" << endl;
  JacobiRecursion jr(a,b,x);
  for (int n(0); n <= 4; n++)
    {
      double exact(0), approx(jr(n));
      switch (n)
	{
	case 0: exact = 1.0; break;
	case 1: exact = 2.0*x; break;
	case 2: exact = 3.75*x*x - 0.75; break;
	case 3: exact = 7.0*x*x*x - 3.0*x; break;
	case 4: exact = 13.125*x*x*x*x -8.75*x*x + 0.625; break;
	}
      cout << "  * T_" << n << "(" << x << ")=" << approx
	   << " (error: " << fabs(approx-exact) << ")" << endl;
    }

  x = 0.1;
  cout << "- Legendre recursion for x=" << x << ":" << endl;
  LegendreRecursion lr(x);
  for (int n(0); n <= 4; n++)
    {
      double exact(0), approx(lr(n));
      switch (n)
	{
	case 0: exact = 1.0; break;
	case 1: exact = x; break;
	case 2: exact = 1.5*x*x - 0.5; break;
	case 3: exact = 2.5*x*x*x - 1.5*x; break;
	case 4: exact = 4.375*x*x*x*x - 3.75*x*x + 0.375; break;
	}
      cout << "  * T_" << n << "(" << x << ")=" << approx
	   << " (error: " << fabs(approx-exact) << ")" << endl;
    }

  Vector<double> coeffs(3);
  coeffs[0] = -18.0;
  coeffs[1] = 4.0;
  coeffs[2] = 1.0; // x^2+4*x-18
  x = 6.0;
  cout << "- forward summation of a polynomial: "
       << MonomialRecursion(x).forwardSummation(coeffs) << endl;
  cout << "- adjoint summation of a polynomial (Horner scheme): "
       << MonomialRecursion(x).adjointSummation(coeffs) << endl;

  return 0;
}

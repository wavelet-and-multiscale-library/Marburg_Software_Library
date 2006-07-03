#include <iostream>
#include <fstream>
#include <numerics/bezier.h>
#include <utils/array1d.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing Bezier stuff from MathTL:: ..." << endl;

  const unsigned int n = 50;

#if 0
  const int d = 3;
  for (int k = -1; k <= d+1; k++) {
    cout << "* some point values of the " << k << "-th Bernstein polynomial of degree " << d << ":" << endl;
    Array1D<double> xi(n+1), yi(n+1);
    for (unsigned int i = 0; i <= n; i++) {
      xi[i] = i/(double)n;
      yi[i] = EvaluateBernsteinPolynomial<d>(k, xi[i]);
    }
    
    cout << "  xi=" << xi << endl;
    cout << "  yi=" << yi << endl;
  }
  for (int k = -1; k <= d+1; k++) {
    cout << "* some point values of the 1st derivative of the " << k << "-th Bernstein polynomial of degree " << d << ":" << endl;
    Array1D<double> xi(n+1), yi(n+1);
    for (unsigned int i = 0; i <= n; i++) {
      xi[i] = i/(double)n;
      yi[i] = EvaluateBernsteinPolynomial_x<d>(k, xi[i]);
    }
    
    cout << "  xi=" << xi << endl;
    cout << "  yi=" << yi << endl;
  }
#endif
  
  double b0 = 0;
  double b1 = 1;
  double b2 = 0;
  double b3 = 0;
  cout << "* some point values from a cubic Bernstein polynomial with Bezier coefficients "
       << "b0=" << b0 << ", b1=" << b1 << ", b2=" << b2 << ", b3=" << b3 << ":" << endl;

  Array1D<double> xi(n+1), yi(n+1);
  for (unsigned int i = 0; i <= n; i++) {
    xi[i] = i/(double)n;
    yi[i] = EvaluateBernsteinPolynomial(b0, b1, b2, b3, xi[i]);
  }

  cout << "  xi=" << xi << endl;
  cout << "  yi=" << yi << endl;
  
  for (unsigned int m = 0; m <= 1; m++) {
    cout << "* some point values from the " << m << "-th cubic Hermite interpolant:" << endl;
    
    for (unsigned int i = 0; i <= n; i++) {
      xi[i] = -1.5 + i*(1.5-(-1.5))/(double)n;
      yi[i] = EvaluateHermiteSpline(m, xi[i]);
    }
    
    cout << "  xi=" << xi << endl;
    cout << "  yi=" << yi << endl;
  }
  
  for (unsigned int m = 0; m <= 1; m++) {
    cout << "* some point values from the first derivative of the " << m << "-th cubic Hermite interpolant:" << endl;
    
    for (unsigned int i = 0; i <= n; i++) {
      xi[i] = -1.5 + i*(1.5-(-1.5))/(double)n;
      yi[i] = EvaluateHermiteSpline_x(m, xi[i]);
    }
    
    cout << "  xi=" << xi << endl;
    cout << "  yi=" << yi << endl;
  }
  
  return 0;
}

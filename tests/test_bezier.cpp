#include <iostream>
#include <fstream>
#include <numerics/bezier.h>

using std::cout;
using std::endl;
using namespace MathTL;

int main()
{
  cout << "Testing Bezier stuff from MathTL:: ..." << endl;

  const int d = 1;
  for (int k = -1; k <= d+1; k++) {
    cout << "* some point values of the " << k << "-th Bernstein polynomial of degree " << d << ":" << endl;
    for (double x = 0.0; x <= 1.0; x+=0.1) {
      cout << "  B_{" << k << "," << d << "}(" << x << ")="
	   <<  EvaluateBernsteinPolynomial<d>(k, x) << endl;
    }
  }
  
  return 0;
}

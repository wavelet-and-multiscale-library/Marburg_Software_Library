#include <iostream>
#include <cmath>
#include <vector>
#include <numerics/ortho_poly.h>
#include <numerics/quadrature.h>
#include <numerics/gauss_quadrature.h>
#include <numerics/b_splines.h>
#include <utils/function.h>
#include <utils/array1d.h>

using std::cout;
using std::endl;

using namespace MathTL;

class TestFunction : public Function<1>
{
public:
  TestFunction() : Function<1>(1) {}
  virtual ~TestFunction(){}
  double value(const Point<1>& p,
	       const unsigned int component = 0) const
  {
    return M_PI*sin(M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p,
		    Vector<double>& values) const
  {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  TestFunction g;

  cout << "Testing GaussRule class ..." << endl;

  LegendrePolynomial leg;
  for (unsigned int N(1); N <= 6; N++)
    {
      GaussRule GL(leg, -1.0, 1.0, N);
      cout << "* Gauss-Legendre rule from the recursion coefficients (N=" << N << "): " << endl;
      Array1D<Point<1> > points;
      Array1D<double> weights;
      GL.get_points(points);
      GL.get_weights(weights);
      cout << "  points: " << points << ", weights: " << weights << endl;
      cout << "  quadrature yields: " << GL.integrate(g) << endl;
    }
 
  return 0;
}

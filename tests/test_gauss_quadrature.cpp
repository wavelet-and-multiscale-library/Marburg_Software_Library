#include <iostream>
#include <cmath>
#include <vector>
#include <numerics/ortho_poly.h>
#include <numerics/quadrature.h>
#include <numerics/gauss_quadrature.h>
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

  for (unsigned int N(1); N <= 8; N++) {
      GaussLegendreRule Gauss(N);
      cout << "* Gauss-Legendre reference rule(N=" << N << "): " << endl;
      Array1D<Point<1> > points;
      Array1D<double> weights;
      Gauss.get_points(points);
      Gauss.get_weights(weights);
      cout << "  points: " << points << ", weights: " << weights << endl;
      cout << "  quadrature yields: " << Gauss.integrate(g) << endl;
    }

  LegendrePolynomial leg;
  for (unsigned int N(1); N <= 10; N++) {
      GaussRule GL(leg, -1.0, 1.0, N);
      cout << "* Gauss-Legendre rule from the recursion coefficients (N=" << N << "): " << endl;
      Array1D<Point<1> > points;
      Array1D<double> weights;
      GL.get_points(points);
      GL.get_weights(weights);
      cout << "  points: " << points << ", weights: " << weights << endl;
      cout << "  quadrature yields: " << GL.integrate(g) << endl;
    }
 
  for (unsigned int N(1); N <= 10; N++) {
#if 1
      // set up 2N+1 monomial moments of w(x)=1 on [-1,1] (yields a bit more exact values)
      Array1D<double> moments(2*N+1);
      for (unsigned int n(0); n < moments.size(); n++)
	if (n%2 == 0)
	  moments[n] = 2./(n+1); // = 1^{n+1}/(n+1) - (-1)^{n+1}/(n+1)
	else
	  moments[n] = 0.0; // odd moments are 0
      
      GaussRule GL(moments, -1.0, 1.0, N);
#else
      // set up 2N+1 monomial moments of w(x)=1 on [0,1]
      Array1D<double> moments(2*N+1);
      for (unsigned int n(0); n < moments.size(); n++)
	moments[n] = 1./(n+1); // = 1^{n+1}/(n+1)
      
      GaussRule GL(moments, 0.0, 1.0, N);
#endif
      cout << "* Gauss-Legendre rule from monomial (N=" << N << "): " << endl;
      Array1D<Point<1> > points;
      Array1D<double> weights;
      GL.get_points(points);
      GL.get_weights(weights);
      cout << "  points: " << points << ", weights: " << weights << endl;
      cout << "  quadrature yields: " << GL.integrate(g) << endl;
    }

  for (unsigned int N(1); N <= 4; N++) {
    // set up generalized moments of w(x)=1 on [-1,1] w.r.t. the (normalized) Legendre polynomials
    Array1D<double> moments(2*N);
    Polynomial<double> w(1);
    for (unsigned int n(0); n < 2*N; n++) {
      moments[n] = (w*leg.assemble(n)).integrate(-1.0, 1.0);
    }
    GaussRule Gleg(moments, leg, -1.0, 1.0, N);

    cout << "* Gauss-Legendre rule from generalized moments"
	 << " (w.r.t. Legendre poly., N="
	 << N << "): " << endl;
    Array1D<Point<1> > points;
    Array1D<double> weights;
    Gleg.get_points(points);
    Gleg.get_weights(weights);
    cout << "  points: " << points << ", weights: " << weights << endl;
    cout << "  quadrature yields: " << Gleg.integrate(g) << endl;

    for (unsigned int power(0); power <= 8; power++)
      {
        Polynomial<double> p;
        p.set_coefficient(power, 1.0);
	double integral = Gleg.integrate(p, Point<1>(-1.0), Point<1>(1.0));
	double exact = (w*p).integrate(-1.0, 1.0, false);
        cout << "  integration of p(x)=" << p << " yields: "
             << integral << ", absolute error: "
             << fabs(integral-exact) << endl;
      }
  }

  return 0;
}

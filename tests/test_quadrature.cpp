#include <iostream>
#include <cmath>
#include <vector>
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

  cout << "Testing QuadratureRule class ..." << endl;
  MidpointRule MPQ;
  cout << "* midpoint rule: " << MPQ.integrate(g) << endl;
//   Array1D<Point<1> > points;
//   Array1D<double> weights;
//   MPQ.get_points(points);
//   MPQ.get_weights(weights);
//   cout << "  (points: " << points << ", weights: " << weights << ")" << endl;
  TrapezoidalRule TQ;
  cout << "* trapezoidal rule: " << TQ.integrate(g) << endl;
  SimpsonRule SQ;
  cout << "* Simpson rule: " << SQ.integrate(g) << endl;
  GaussLegendreRule<1> G1Q;
  cout << "* Gauss-Legendre rule (N=1): "
       << G1Q.integrate(g) << endl;
  GaussLegendreRule<2> G2Q;
  cout << "* Gauss-Legendre rule (N=2): "
       << G2Q.integrate(g) << endl;
  GaussLegendreRule<3> G3Q;
  cout << "* Gauss-Legendre rule (N=3): "
       << G3Q.integrate(g) << endl;
  GaussLegendreRule<4> G4Q;
  cout << "* Gauss-Legendre rule (N=4): "
       << G4Q.integrate(g) << endl;
  GaussLegendreRule<5> G5Q;
  cout << "* Gauss-Legendre rule (N=5): "
       << G5Q.integrate(g) << endl;
  GaussLegendreRule<6> G6Q;
  cout << "* Gauss-Legendre rule (N=6): "
       << G6Q.integrate(g) << endl;
  ClosedNewtonCotesRule<1> CNC1Q;
  cout << "* closed Newton-Cotes rule (N=1): "
       << CNC1Q.integrate(g) << endl;
  ClosedNewtonCotesRule<2> CNC2Q;
  cout << "* closed Newton-Cotes rule (N=2): "
       << CNC2Q.integrate(g) << endl;
  ClosedNewtonCotesRule<3> CNC3Q;
  cout << "* closed Newton-Cotes rule (N=3): "
       << CNC3Q.integrate(g) << endl;
  ClosedNewtonCotesRule<4> CNC4Q;
  cout << "* closed Newton-Cotes rule (N=4): "
       << CNC4Q.integrate(g) << endl;
  ClosedNewtonCotesRule<5> CNC5Q;
  cout << "* closed Newton-Cotes rule (N=5): "
       << CNC5Q.integrate(g) << endl;
  ClosedNewtonCotesRule<10> CNC10Q;
  cout << "* closed Newton-Cotes rule (N=10): "
       << CNC10Q.integrate(g) << endl;
  ClosedNewtonCotesRule<20> CNC20Q;
  cout << "* closed Newton-Cotes rule (N=20): "
       << CNC20Q.integrate(g) << endl;

  cout << "- check integration over arbitrary intervals, e.g. [0,3]: " << endl;
  cout << "* Simpson rule: "
       << SQ.integrate(g, Point<1>(0.0), Point<1>(3.0)) << endl;
  cout << "* Gauss-Legendre rule (N=1): "
       << G1Q.integrate(g, Point<1>(0.0), Point<1>(3.0)) << endl;
  cout << "* Gauss-Legendre rule (N=2): "
       << G2Q.integrate(g, Point<1>(0.0), Point<1>(3.0)) << endl;
  cout << "* Gauss-Legendre rule (N=3): "
       << G3Q.integrate(g, Point<1>(0.0), Point<1>(3.0)) << endl;
  cout << "* Gauss-Legendre rule (N=4): "
       << G4Q.integrate(g, Point<1>(0.0), Point<1>(3.0)) << endl;
  cout << "* Gauss-Legendre rule (N=5): "
       << G5Q.integrate(g, Point<1>(0.0), Point<1>(3.0)) << endl;
  cout << "* Gauss-Legendre rule (N=6): "
       << G6Q.integrate(g, Point<1>(0.0), Point<1>(3.0)) << endl;
  
  QuadratureRule<2> TensorQ(SQ, G2Q);
  Array1D<Point<2> > points2D;
  Array1D<double> weights2D;
  TensorQ.get_points(points2D);
  TensorQ.get_weights(weights2D);
  cout << "- construct a 2D Quadrature rule as a tensor product," << endl
       << "  points: " << points2D << ", weights: " << weights2D << endl;
  
  cout << "- construct 1D composite rules:" << endl;
  for (unsigned int N(1); N <= (1<<4); N <<= 1)
    {
      CompositeRule<1> CompTQ(TQ, N);
      cout << "* trapezoidal composite rule with " << N
	   << " subdivisions," << endl
	   << "  integration: "
	   << CompTQ.integrate(g) << endl;
//       Array1D<Point<1> > points;
//       Array1D<double> weights;
//       CompTQ.get_points(points);
//       CompTQ.get_weights(weights);
//       cout << "  points: " << points << ", weights: " << weights << endl;
    }

  for (unsigned int N(1); N <= (1<<4); N <<= 1)
    {
      CompositeRule<1> CompG2Q(G2Q, N);
      cout << "* Gauss 2-point composite rule with " << N
 	   << " subdivisions," << endl
	   << "  integration: "
	   << CompG2Q.integrate(g) << endl;
//       Array1D<Point<1> > points;
//       Array1D<double> weights;
//       CompG2Q.get_points(points);
//       CompG2Q.get_weights(weights);
//       cout << "  points: " << points << ", weights: " << weights << endl;
    }

  Bspline<2> N2;
  cout << "- integrating N_2 with composite rules:" << endl;
  for (int N(1); N <= (1<<3); N <<= 1)
    {
      cout << "* composite rule (trapez., N=" << N << "): "
	   << CompositeRule<1>(TQ, N).integrate(N2, Point<1>(0.0), Point<1>(2.0))
 	   << endl;
    }
  
  return 0;
}

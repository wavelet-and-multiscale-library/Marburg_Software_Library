#include <iostream>
#include <cmath>
#include <vector>
#include <numerics/quadrature.h>
#include <numerics/gauss_quadrature.h>
#include <numerics/cardinal_splines.h>
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

class TestFunction_stretched : public Function<1>
{
public:
  TestFunction_stretched() : Function<1>(1) {}
  virtual ~TestFunction_stretched(){}
  double value(const Point<1>& p,
	       const unsigned int component = 0) const
  {
    return M_PI*sin(M_PI*p[0]/3.)/3.;
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
  TestFunction_stretched g_stretched;

  cout << "Testing QuadratureRule class ..." << endl;

  MidpointRule MPQ;
  cout << "* midpoint rule: " << MPQ.integrate(g) << endl;

  TrapezoidalRule TQ;
  cout << "* trapezoidal rule: " << TQ.integrate(g) << endl;

  SimpsonRule SQ;
  cout << "* Simpson rule: " << SQ.integrate(g) << endl;

  for (unsigned int N(1); N <= 6; N++)
    {
      GaussLegendreRule Gauss(N);
      cout << "* Gauss-Legendre rule (N=" << N << "): "
	   << Gauss.integrate(g) << endl;
    }
  
  for (unsigned int N(1); N <= 5; N++)
    {
      ClosedNewtonCotesRule NewtonCotes(N);
      cout << "* closed Newton-Cotes rule (N=" << N << "): "
	   << NewtonCotes.integrate(g) << endl;
    }
  for (unsigned int N(10); N <= 25; N+=5)
    {
      ClosedNewtonCotesRule NewtonCotes(N);
      cout << "* closed Newton-Cotes rule (N=" << N << "): "
	   << NewtonCotes.integrate(g) << endl;
    }

  cout << "- check integration over arbitrary intervals, e.g. [0,3]: " << endl;
  cout << "* Simpson rule: "
       << SQ.integrate(g_stretched, Point<1>(0.0), Point<1>(3.0)) << endl;

  for (unsigned int N(1); N <= 6; N++)
    {
      GaussLegendreRule Gauss(N);
      cout << "* Gauss-Legendre rule (N=" << N << "): "
	   << Gauss.integrate(g_stretched, Point<1>(0.0), Point<1>(3.0)) << endl;
    }
  
  QuadratureRule<2> TensorQ(SQ, GaussLegendreRule(2));
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
      CompositeRule<1> CompG2Q(GaussLegendreRule(2), N);
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

  CardinalBSpline<2> N2;
  cout << "- integrating N_2 with composite rules (should be 1):" << endl;
  for (int N(1); N <= (1<<3); N <<= 1)
    {
      cout << "* composite rule (trapez., N=" << N << "): "
	   << CompositeRule<1>(TQ, N).integrate(N2, Point<1>(0.0), Point<1>(2.0))
 	   << endl;
    }
  
  CardinalBSpline<3> N3;
  cout << "- integrating N_3 with composite rules (should be 1):" << endl;
  for (int N(1); N <= (1<<5); N <<= 1)
    {
      cout << "* composite rule (trapez., N=" << N << "): "
	   << CompositeRule<1>(TQ, N).integrate(N3, Point<1>(0.0), Point<1>(3.0))
 	   << endl;
    }
  
  CardinalBSpline<4> N4;
  cout << "- integrating N_4 with composite rules (should be 1):" << endl;
  for (int N(1); N <= (1<<5); N <<= 1)
    {
      cout << "* composite rule (trapez., N=" << N << "): "
	   << CompositeRule<1>(TQ, N).integrate(N4, Point<1>(0.0), Point<1>(4.0))
 	   << endl;
    }
  
  return 0;
}

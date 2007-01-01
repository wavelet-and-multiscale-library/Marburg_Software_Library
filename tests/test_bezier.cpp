#include <iostream>
#include <fstream>
#include <geometry/sampled_mapping.h>
#include <numerics/bezier.h>
#include <numerics/quadrature.h>
#include <numerics/gauss_quadrature.h>
#include <utils/array1d.h>
#include <utils/function.h>

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

  Array1D<double> xi(n+1), yi(n+1), zi(n+1), err(n+1);
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
      
      if (xi[i] <= -1 || xi[i] >= 1)
	zi[i] = 0;
      else {
	if (xi[i] <= 0) {
	  if (m == 0)
	    zi[i] = (xi[i]+1)*(xi[i]+1)*(1-2*xi[i]);
	  else
 	    zi[i] = xi[i]*(xi[i]+1)*(xi[i]+1);
	} else {
	  zi[i] = 0;
 	  if (m == 0)
	    zi[i] = (1-xi[i])*(1-xi[i])*(2*xi[i]+1);
 	  else
 	    zi[i] = xi[i]*(xi[i]-1)*(xi[i]-1);
	}
      }

      err[i] = fabs(yi[i]-zi[i]);
    }
    
    cout << "  xi=" << xi << endl;
    cout << "  yi=" << yi << endl;
    cout << "  zi=" << zi << endl;
    cout << "  error=" << err << endl;
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

  cout << "- examine the two components of a Hermite interpolatory spline:" << endl;
  CubicHermiteInterpolant_td phi1(0, 0, 0), phi2(0, 0, 1);
  CompositeRule<1> qrule(GaussLegendreRule(4), 4);
  ProductFunction<1> integrand(&phi1, &phi1);
  cout << "||phi_0||^2_2=" << qrule.integrate(integrand, Point<1>(-1), Point<1>(1))
       << " (expected " << 26./35. << ")" << endl;
  ProductFunction<1> integrand2(&phi1, &phi2);
  cout << "<phi_0,phi_0(.-1)>=" << qrule.integrate(integrand2, Point<1>(-1), Point<1>(1))
       << " (expected " << 9./70. << ")" << endl;
  phi2.set_c(1); phi2.set_k(0);
  ProductFunction<1> integrand3(&phi1, &phi2);
  cout << "<phi_0,phi_1>=" << qrule.integrate(integrand3, Point<1>(-1), Point<1>(1))
       << " (expected " << 0. << ")" << endl;
  ProductFunction<1> integrand4(&phi2, &phi2);
  cout << "||phi_1||^2_2=" << qrule.integrate(integrand4, Point<1>(-1), Point<1>(1))
       << " (expected " << 2./105. << ")" << endl;
  phi2.set_k(1);
  ProductFunction<1> integrand5(&phi1, &phi2);
  cout << "<phi_0,phi_1(.-1)>=" << qrule.integrate(integrand5, Point<1>(-1), Point<1>(1))
       << " (expected " << -13./420. << ")" << endl;
  phi1.set_c(1);
  ProductFunction<1> integrand6(&phi1, &phi2);
  cout << "<phi_1,phi_1(.-1)>=" << qrule.integrate(integrand6, Point<1>(-1), Point<1>(1))
       << " (expected " << -1./140. << ")" << endl;
  
#if 1
  cout << "* plot a cubic Hermite interpolant into a file..." << endl;
  CubicHermiteInterpolant2D_td fh(3, 1, 1, 0, 0);
  Grid<2> fh_grid(Point<2>(-1,0), Point<2>(0,1), 1<<5);
  SampledMapping<2> fh_sm(fh_grid, fh);
  std::ofstream fh_stream("hermite_spline.m");
  fh_sm.octave_output(fh_stream);
  fh_stream.close();
#endif

  return 0;
}

#include <iostream>
#include <fstream>
#include <sstream>
#include <algebra/matrix.h>
#include <geometry/point.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>
#include <numerics/cardinal_splines.h>
#include <numerics/splines.h>
#include <numerics/schoenberg_splines.h>
#include <utils/array1d.h>

using namespace std;
using namespace MathTL;

int main()
{
  cout << "Testing MathTL::(Cardinal)Spline ..." << endl;

  cout << "* some point values of the linear B-spline:" << endl;
  for (double x = -0.5; x <= 2.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  EvaluateCardinalBSpline<2>(0, x) << endl;
  }

  cout << "* some point values of a shifted quadratic B-spline:" << endl;
  for (double x = -0.5; x <= 3.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  EvaluateCardinalBSpline(3, 1, x) << endl;
  }  

  cout << "* some point values of the CDF function phi, d=2:" << endl;
  for (double x = -1.5; x <= 2.5; x+=0.5) {
    cout << "  phi(" << x << ")="
	 <<  EvaluateCardinalBSpline_td(2, 0, 0, x) << endl;
  }  
  
  cout << "* some point values of phi_{j,k}, d=2, j=1, k=0:" << endl;
  for (double x = -1.5; x <= 2.5; x+=0.5) {
    cout << "  phi(" << x << ")="
	 <<  EvaluateCardinalBSpline_td(2, 1, 0, x) << endl;
  }  

  cout << "- some point values of phi_{3,4}, d=2:" << endl;
  for (double x = 0.0; x <= 1.0; x += ldexp(1.0, -5))
    {
      cout << "phi(" << x << ")="
	   << EvaluateCardinalBSpline_td(2, 3, 4, x) << endl;
    }
  
  cout << "* some point values of the default constant spline:" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  Spline<1>().value(Point<1>(x)) << endl;
  }

  cout << "* some point values of the default linear spline:" << endl;
  for (double x = -0.5; x <= 2.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  Spline<2>().value(Point<1>(x)) << endl;
  }

  cout << "* some point values of the default quadratic spline:" << endl;
  for (double x = -0.5; x <= 3.5; x+=0.5) {
    cout << "  N(" << x << ")="
	 <<  Spline<3>().value(Point<1>(x)) << endl;
  }

#if 0
  cout << "* a linear spline from a given knot sequence:" << endl;
  Array1D<double> knots(4+2);
  knots[0] =  0.0;
  knots[1] =  1.0;
  knots[2] =  1.5;
  knots[3] =  2.0;
  knots[4] =  2.5;
  knots[5] =  3.0;
  Array1D<double> coeffs(4);
  coeffs[0] = 1;
  coeffs[1] = 0;
  coeffs[2] = 0;
  coeffs[3] = 0;
  SampledMapping<1> sm1(Grid<1>(0.0, 3.0, 30), Spline<2>(knots, coeffs));
  sm1.matlab_output(cout);
#endif

#if 0
  cout << "* writing point values of a linear spline with multiple knots to a file..." << endl;
  knots[0] =  0.0;
  knots[1] =  0.0;
  knots[2] =  1.0;
  knots[3] =  2.0;
  knots[4] =  3.0;
  knots[5] =  3.0;
  coeffs[0] = 1;
  coeffs[1] = 0;
  coeffs[2] = 1;
  coeffs[3] = 0;
  std::ofstream fs("multipleknots2.m");
  SampledMapping<1> sm2(Grid<1>(0.0, 3.0, 60), Spline<2>(knots, coeffs));
  sm2.matlab_output(fs);
  fs.close();
#endif

#if 0
  {
    const unsigned int d = 3;
    Array1D<double> knots(5+d);
    cout << "* writing point values of various quadratic splines with multiple knots to a file..." << endl;
    knots[0] =  0.0;
    knots[1] =  0.0;
    knots[2] =  1.0;
    knots[3] =  1.0;
    knots[4] =  2.0;
    knots[5] =  2.0;
    knots[6] =  3.0;
    knots[7] =  3.0;
    Array1D<double> coeffs(5);
    for (unsigned int i = 0; i < coeffs.size(); i++) {
      for (unsigned int j = 0; j < coeffs.size(); j++) {
	coeffs[j] = (i == j ? 1.0 : 0.0);
      }
      ostringstream filename;
      filename << "spline_" << i << ".m";
      std::ofstream fs(filename.str().c_str());
      SampledMapping<1> sm(Grid<1>(0.0, 3.0, 60), Spline<d>(knots, coeffs));
      sm.matlab_output(fs);
      fs.close();
    }
  }
#endif

#if 0
  cout << "* writing point values of a cubic spline with multiple knots to a file..." << endl;
  knots.resize(7);
  knots[0] =  0.0;
  knots[1] =  0.0;
  knots[2] =  0.0;
  knots[3] =  1.0;
  knots[4] =  2.0;
  knots[5] =  2.0;
  knots[6] =  2.0;
  coeffs.resize(4);
  coeffs[0] = 1;
  coeffs[1] = 0;
  coeffs[2] = 1;
  coeffs[3] = 0;
  std::ofstream fs2("multipleknots3.m");
  SampledMapping<1> sm3(Grid<1>(0.0, 2.0, 1000), Spline<3>(knots, coeffs));
  sm3.matlab_output(fs2);
  fs2.close();
#endif

  cout << "- testing some Schoenberg splines:" << endl;
  cout << "* some point values of the 0-th first-order spline:" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.1) {
    cout << "  N_{0,1}(" << x << ")="
	 <<  EvaluateSchoenbergBSpline<1>(0,x) << endl;
  }
  cout << "* some point values of the (-1)-th second-order spline:" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.1) {
    cout << "  N_{-1,2}(" << x << ")="
	 <<  EvaluateSchoenbergBSpline<2>(-1,x) << endl;
  }

  SchoenbergBSpline_td<2> sbs(0, -1);
  cout << "* compare with the values from a function object:" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.1) {
    cout << "  N_{-1,2}(" << x << ")="
	 <<  sbs.value(Point<1>(x)) << endl;
  }
  
  cout << "* some point values of its derivative:" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.1) {
    cout << "  N_{-1,2}'(" << x << ")="
	 <<  EvaluateSchoenbergBSpline_x<2>(-1,x) << endl;
  }
  cout << "* some point values of the derivative of N_{0,2}:" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.1) {
    cout << "  N_{0,2}'(" << x << ")="
	 <<  EvaluateSchoenbergBSpline_x<2>(0,x) << endl;
  }

  const int d = 3;

  CardinalKnotSequence<d> cknots;
  SchoenbergKnotSequence<d> sknots;
  cout << "- (begin of) the cardinal knot sequence for d=" << d << ":" << endl;
  for (int k = cknots.k0(); k <= 3; k++)
    cout << "  t_{" << k << "}=" << cknots.knot(k) << endl;
  cout << "- (begin of) the Schoenberg knot sequence for d=" << d << ":" << endl;
  for (int k = sknots.k0(); k <= 3; k++)
    cout << "  t_{" << k << "}=" << sknots.knot(k) << endl;

  cout << "- some point values of the (-1)-th Schoenberg spline of order d=" << d << ":" << endl;
  for (double x = -0.5; x <= 1.5; x+=0.1) {
    cout << "  N_{-1," << d << "}(" << x << ")="
	 <<  evaluate_Bspline<d>(&sknots,-1,x) << endl;
  }
  
  Matrix<double> M;
  compute_Bspline_refinement_matrix<d>(&cknots, M);
  cout << "- refinement matrix for the cardinal boundary B-splines (d=" << d << "):" << endl
       << M;
  compute_Bspline_refinement_matrix<d>(&sknots, M);
  cout << "- refinement matrix for the Schoenberg boundary B-splines (d=" << d << "):" << endl
       << M;
}

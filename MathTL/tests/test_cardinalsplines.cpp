#include <iostream>
#include <fstream>
#include <sstream>
#include <algebra/matrix.h>
#include <algebra/piecewise.h>
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






  cout << "* some point values of phi_{j,k}, d=2, j=1, k=0:" << endl;
  for (double x = -1.5; x <= 2.5; x+=0.2) {
    cout << "  phi(" << x << ")="
	 <<  EvaluateCardinalBSpline_td(2, 0, 1.5, x) << endl;
  } 



  cout << "* some point values of phi_{j,k}, d=2, j=1, k=0: with expand as picewise" << endl;
  for (double x = -1.5; x <= 2.5; x+=0.2) {
    Piecewise<double> p;
    p = ExpandBspline<3>(1, -1);
    cout << "  phi(" << x << ")="
	 <<  p(x) << endl;
  } 



  cout << "* some point values of Schoenberg Primal [P] function phi_{j,k}(x) = 2^{j/2}N_{k-d/2,d}(2^jx):" << endl;
  for (double x = -0.6; x <= 1.6; x+=0.2){
    cout << "  N_{-1,2}'(" << x << ")="
	 <<  EvaluateSchoenbergBSpline_td_xx<2>(2, 5,x) << endl;
  }


cout << "* some point values of Schoenberg Primal [P] function phi_{j,k}(x) = 2^{j/2}N_{k-d/2,d}(2^jx) with expand as picewise:" << endl;
  for (double x = -0.6; x <= 1.6; x+=0.2) {
    Piecewise<double> p;
    p = ExpandSchoenbergBspline<2>(2, 5,0);
    cout << "  N(" << x << ")="
	 <<  p.secondDerivative(x) << endl;
  }

#if 0
  cout << "* some point values of Schoenberg Primal [P] function phi_{j,k}(x) = 2^{j/2}N_{k-d/2,d}(2^jx):" << endl;
    cout <<  EvaluateSchoenbergBSpline_td<4>(0, -1, 0.0000) << endl;
    cout <<  EvaluateSchoenbergBSpline_td_x<4>(0, -1, 0.0000) << endl;
    cout <<  EvaluateSchoenbergBSpline_td_xx<4>(0, -1, 0.0000) << endl;
    double j,a,b,c,d;    
    j = (EvaluateSchoenbergBSpline_td_xx<4>(0, -1, 0.00002) - EvaluateSchoenbergBSpline_td_xx<4>(0, -1, 0.00001))/(0.00001);
    cout <<  j << endl;
	a= 1<<5;
	cout << "a=" << a  << endl;	
  



  //double EvaluateSchoenbergBSpline_td(const int j, const int k, const double x)
  Array1D<double> values1(80);
  for (double x = -1.0; x <= 3; x+=0.05) {
	values1[(int)floor((x+1.0)*20)] = EvaluateSchoenbergBSpline_td<2>(0, 1, x);
  }
    std::ofstream fs1("Schoenberg1(0,0.5).m");
    SampledMapping<1> sm1(Grid<1>(-1.0, 3.0, 80), values1);
    sm1.matlab_output(fs1);
    fs1.close();


  //EvaluateCardinalBSpline_td(const int j, const int k, const double x)
  Array1D<double> values2(80);
  for (double x = -1.0; x <= 3; x+=0.05) {
	values2[(int)floor((x+1.0)*20)] = EvaluateCardinalBSpline_td<2>(0, 1, x);
  }
    std::ofstream fs2("Cardinal1(0,0.5).m");
    SampledMapping<1> sm2(Grid<1>(-1.0, 3.0, 80), values2);
    sm2.matlab_output(fs2);
    fs2.close();

#endif


}

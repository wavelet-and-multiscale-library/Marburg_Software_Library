#include <iostream>
#include <fstream>
#include <cmath>
#include <numerics/goertzel_reinsch.h>
#include <algebra/vector.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

using std::cout;
using std::endl;

using namespace MathTL;

int main()
{
  cout << "Testing the Goertzel-Reinsch algorithm..." << endl;

#if 0
  double x;

  Vector<double> coeffs(1001);
  coeffs[1000] = 1;
  
  {
    x = M_PI/2000.;
    
    Goertzel sg(coeffs, true);
    Goertzel cg(coeffs, false);
    
    GoertzelReinsch sr(coeffs, true);
    GoertzelReinsch cr(coeffs, false);
    
    const double exacts = sin(1000*x);
    cout << "sin(1000*" << x << ")=" << exacts << endl;
    cout << "s_G(" << x << ")=" << sg.value(Point<1>(x)) << ", "
 	 << "abs. error " << fabs(sg.value(Point<1>(x))-exacts) << endl;
    cout << "s_R(" << x << ")=" << sr.value(Point<1>(x)) << ", "
 	 << "abs. error " << fabs(sr.value(Point<1>(x))-exacts) << endl;
    
    const double exactc = cos(1000*x);
    cout << "cos(1000*" << x << ")=" << exactc << endl;
    cout << "c_G(" << x << ")=" << cg.value(Point<1>(x)) << ", "
	 << "abs. error " << fabs(cg.value(Point<1>(x))-exactc) << endl;
    cout << "c_R(" << x << ")=" << cr.value(Point<1>(x)) << ", "
	 << "abs. error " << fabs(cr.value(Point<1>(x))-exactc) << endl;
  }
#endif

//   typedef float C; // leads to pointwise difference of approx. 1e-4
  typedef double C; // no pointwise difference visible

  const unsigned int N = 100;
  Vector<C> fewcoeffs(N);
  fewcoeffs[N-1] = -1.0;
  
  cout << "a coefficient set: " << fewcoeffs << endl;
  cout << "plotting the corresponding sine series to a file..." << endl;
  Goertzel<C> sineseries_g(fewcoeffs);
  std::ofstream fs_g("sineseries_g.m");
  SampledMapping<1,C> sm_g(Grid<1>(0.0, M_PI, 10000), sineseries_g);
  sm_g.matlab_output(fs_g);
  fs_g.close();
  GoertzelReinsch<C> sineseries_r(fewcoeffs);
  std::ofstream fs_r("sineseries_r.m");
  SampledMapping<1,C> sm_r(Grid<1>(0.0, M_PI, 10000), sineseries_r);
  sm_r.matlab_output(fs_r);
  fs_r.close();
  cout << "... done!" << endl;

  return 0;
}

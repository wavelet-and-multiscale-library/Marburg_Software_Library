#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algebra/infinite_vector.h>
#include <algebra/matrix.h>
#include <geometry/grid.h>
#include <interval/periodic.h>
#include <utils/function.h>
#include <ring/ring_basis.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

// the constant function f(x)=23
class RingFunction1
  : public Function<2>
{
public:
  inline double value(const Point<2>& p, const unsigned int component = 0) const {
    return 23;
  }
  void vector_value(const Point<2> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

// the radial hat function
class RingFunction2
  : public Function<2>
{
public:
  inline double value(const Point<2>& p, const unsigned int component = 0) const {
    const double r = sqrt(p[0]*p[0]+p[1]*p[1]);
    return max(1-fabs(1.5-r),0.);
  }
  void vector_value(const Point<2> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};

int main()
{
  cout << "Testing wavelet bases on the ring-shaped domain ..." << endl;

  const int d  = 2;
  const int dt = 2;

  typedef RingBasis<d,dt,0,0> Basis;
  Basis basis;

  typedef Basis::Index Index;

  Function<2>* u = 0;
  const int testcase = 2;
  switch(testcase) {
  case 2:
    u = new RingFunction2();
    break;
  default: // also for testcase 1
    u = new RingFunction1();
  }

#if 0
  // plot the test function
  cout << "- plotting the test function..." << endl;
  const int N = 50;
  Matrix<double> gridx(N+1), gridy(N+1);
  for (int i = 0; i <= N; i++) {
    const double phi = i/(double)N;
    for (int j = 0; j <= N; j++) {
      const double r = 1.0+j/(double)N;
      gridx.set_entry(j, i, r*cos(2*M_PI*phi));
      gridy.set_entry(j, i, r*sin(2*M_PI*phi));
    }
  }
  Grid<2> grid(gridx, gridy);
  SampledMapping<2> rf(grid, *u);
  std::ofstream rfstream("test_function.m");
  rf.matlab_output(rfstream);
  rfstream.close();
  cout << "  ...done, see file test_function.m!" << endl;
#endif
  
  InfiniteVector<double,Index> coeffs;
  
  const int j0 = basis.j0();
  const int jmax = j0+1;
  
#if 0
  basis.expand(u, true, j0, coeffs);
  cout << "- integrals of the test function u against all wavelets up to level jmax="
       << j0 << ":" << endl
       << coeffs << endl;
#endif

#if 1
  basis.expand(u, false, j0, coeffs);
  cout << "- approximate expansion coefficients of the test function in the primal basis:" << endl
       << coeffs << endl;
#endif
  
  
  if (u) delete u;
  
  return 0;
}

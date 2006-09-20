#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <list>

#include <time.h>

#define _MATHTL_ONESTEPSCHEME_VERBOSITY 1
#define _WAVELETTL_CDD1_VERBOSITY 0

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>
#include <geometry/grid.h>
#include <geometry/sampled_mapping.h>

#include <interval/i_index.h>
#include <interval/ds_basis.h>
#include <interval/ds_expansion.h>
#include <interval/p_basis.h>
#include <interval/p_expansion.h>
#include <galerkin/sturm_equation.h>
#include <galerkin/cached_problem.h>

#include <parabolic/row_stage_equation.h>

using namespace std;
using namespace MathTL;
using namespace WaveletTL;

// 1d Poisson equation -u''=g, u(0)=u(1)=0
class Poisson1D
  : public SimpleSturmBVP
{
public:
  double p(const double t) const {
    return 1;
  }
  double p_prime(const double t) const {
    return 0;
  }
  double q(const double t) const {
    return 0;
  }
  double g(const double t) const {
    return 0; // dummy rhs
  }
  bool bc_left() const { return true; }
  bool bc_right() const { return true; }
};

class eigenfunction : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return M_PI*M_PI*sin(M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};


int main()
{
  cout << "Testing ROWStageEquation..." << endl;

  // setup elliptic operator -Delta
  Poisson1D poisson_bvp;

  // setup 1D wavelet basis
  const int d  = 3;
  const int dT = 3;
//   typedef DSBasis<d,dT> Basis;
  typedef PBasis<d,dT> Basis;

  const int jmax = 5;

  typedef Basis::Index Index;
  typedef InfiniteVector<double,Index> V;

  // setup elliptic operator equation, reformulated in ell_2
  typedef SturmEquation<Basis> EllipticEquation;
  EllipticEquation poisson_equation(poisson_bvp, false);

  // select a ROW method
  ROWMethod<Vector<double> > row_method(WMethod<Vector<double> >::ROS2);

  // setup the ROW stage equation object
  ROWStageEquation<EllipticEquation> row_stage_equation
    (&row_method, &poisson_equation, 0, 0, jmax);

  //
  //
  // handle several test cases:
  // 1: u0 = 0, f(t,x) = pi^2*sin(pi*x), u(t,x)=(1-exp(-pi^2*t))*sin(pi*x)

#define _TESTCASE 1

  //
  //
  // setup the initial value u0
  V u0;

#if _TESTCASE == 1
  // do nothing, u0=0
#endif

  //
  //
  // setup driving term f

#if _TESTCASE == 1
  Function<1>* f = new eigenfunction();
#endif

  //
  //
  // setup temporal derivative of driving term f

#if _TESTCASE == 1
  Function<1>* ft = 0;
#endif




  
  if (f) delete f;
  if (ft) delete ft;
  
  return 0;
}

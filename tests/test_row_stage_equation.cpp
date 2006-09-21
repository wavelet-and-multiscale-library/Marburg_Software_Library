#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <list>

#include <time.h>

#define _MATHTL_ONESTEPSCHEME_VERBOSITY 1
#define _WAVELETTL_CDD1_VERBOSITY 0

#include <algebra/infinite_vector.h>
#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/matrix.h>
#include <numerics/sturm_bvp.h>
#include <numerics/iteratsolv.h>
#include <numerics/w_method.h>
#include <numerics/row_method.h>
#include <numerics/one_step_scheme.h>
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
#include <adaptive/cdd1.h>

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

class exact_solution_1 : public Function<1> {
public:
  inline double value(const Point<1>& p, const unsigned int component = 0) const {
    return (1-exp(-M_PI*M_PI*get_time()))*sin(M_PI*p[0]);
  }
  
  void vector_value(const Point<1> &p, Vector<double>& values) const {
    values.resize(1, false);
    values[0] = value(p);
  }
};



/*!
  given a ROW method, t_n and an approximation Du^{(n)},
  compute an approximation of Du^{(n+1)} for the stepwidth h;
  compute also an error estimate
*/
template <class ELLIPTIC>
void increment(const ELLIPTIC* elliptic,
	       ROWStageEquation<ELLIPTIC>* row_stage_equation,
	       const double tolerance,
	       const double t_n,
	       const double h,
	       const InfiniteVector<double,typename ELLIPTIC::Index>& D_un,
	       InfiniteVector<double,typename ELLIPTIC::Index>& D_unplus1,
	       InfiniteVector<double,typename ELLIPTIC::Index>& error_estimate,
	       const int jmax)
{
  typedef typename ELLIPTIC::Index Index;
  typedef InfiniteVector<double,Index> V; // without "typename"...

  const unsigned int stages = row_stage_equation->row_method_->A.row_dimension(); // for readability

  std::list<V> Dalpha_uis; // will hold D_alpha u_i, i=1,...,s

  // solve the stage equations sequentially

  const double alpha = 1./(h*row_stage_equation->row_method_->C(0,0)); // we assume that the gamma_{i,i} are all the same
  cout << "increment: alpha=" << alpha << endl;
  row_stage_equation->set_alpha(alpha);

  for (unsigned int i(0); i < stages; i++) {  
    // setup right-hand side of i-th stage equation
    cout << "increment: setup rhs of stage #" << i << ":" << endl;
    row_stage_equation->setup_rhs(i, tolerance, t_n, h, D_un, Dalpha_uis, jmax);
    cout << "  ... done" << endl;

    // solve i-th stage equation
    V result;
    cout << "increment: call CDD1_SOLVE()..." << endl;
    CDD1_SOLVE(*row_stage_equation, tolerance,
	       D_un, // guess
	       result,
	       1.0,
	       1.0,
	       jmax); // B_alpha*D_alpha x = D_alpha^{-1}y
    cout << "increment: CDD1_SOLVE() done" << endl;

    Dalpha_uis.push_back(result);
  }

  // aggregate the stage solutions for an approximation of Du^{(n+1)}
  D_unplus1.clear();
  error_estimate.clear();
  typename std::list<V>::const_iterator it = Dalpha_uis.begin();
  for (unsigned int i(0); i < stages; i++, ++it) {
    D_unplus1.add(row_stage_equation->row_method_->m[i], *it);
    error_estimate.add(row_stage_equation->row_method_->e[i], *it);
  }
  D_unplus1.scale(row_stage_equation, -1);      // D_unplus1 *= D_alpha^{-1}
  D_unplus1.scale(elliptic, 1);                 // D_unplus1 *= D
  D_unplus1.add(D_un);
  error_estimate.scale(row_stage_equation, -1); // dito
  error_estimate.scale(elliptic, 1);            //
}

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

  const int jmax = 6;

  typedef Basis::Index Index;
  typedef InfiniteVector<double,Index> V;

  // setup elliptic operator equation, reformulated in ell_2
  typedef SturmEquation<Basis> EllipticEquation;
  EllipticEquation poisson_equation(poisson_bvp, false); // do not compute the rhs
//   CachedProblem<SturmEquation<Basis> > celliptic(&poisson_equation,  2.3612 , 13.3116); // PBasis d=2,dT=2
  CachedProblem<SturmEquation<Basis> > celliptic(&poisson_equation,  1.87567, 6.78415); // PBasis d=3,dT=3

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

  //
  //
  // setup exact solution (if available)

#if _TESTCASE == 1
  exact_solution_1 uexact;
#endif

  //
  //
  // select a ROW method
  ROWMethod<Vector<double> > row_method(WMethod<Vector<double> >::ROS2);

  //
  //
  // setup the ROW stage equation object
  ROWStageEquation<CachedProblem<EllipticEquation> > row_stage_equation(&row_method, &celliptic, f, ft);


#if 1
  // nonadaptive solution with constant temporal stepsize h

  V temp, error_estimate, result;
  IVPSolution<V> results;

  const int expo1 = 4;
  const int expo2 = expo1;

  for (int expo = expo1; expo <= expo2; expo++) {
    temp = u0;
    temp.scale(&celliptic, 1); // temp *= D
    
    const int resolution = 10;
    SampledMapping<1> u0_plot(evaluate(celliptic.basis(), u0, true, resolution));

    std::ofstream resultstream;
    resultstream.open("u0.m");
    u0_plot.matlab_output(resultstream);
    resultstream.close();

    int N = 1<<expo;
    double h = 1.0/N;
    cout << "  h = " << h << ":" << endl;

    results.t.push_back(0);
    results.u.push_back(temp);

    for (int i = 1; i <= N; i++) {
      cout << "---------------- before increment " << i << " -----------------------" << endl;

      increment(&celliptic, &row_stage_equation, 1e-5, (i-1)*h, h, temp, result, error_estimate, jmax);
      results.t.push_back(i*h);
      results.u.push_back(result);
      temp = result;

      ostringstream output_filename;
      output_filename << "u" << i << ".m";
      resultstream.open(output_filename.str().c_str());
      result.scale(&celliptic, -1);
      SampledMapping<1> ui_plot(evaluate(celliptic.basis(), result, true, resolution));
      result.scale(&celliptic, 1);
      ui_plot.matlab_output(resultstream);
      resultstream.close();

      cout << "---------------- after increment() -----------------------" << endl;

      V uexact_coeffs;
      uexact.set_time(i*h);
      expand(&uexact, celliptic.basis(), false, jmax, uexact_coeffs);
      result.scale(&celliptic, -1);
      cout << "  ell_2 error at t=" << i*h << ": " << l2_norm(result - uexact_coeffs) << endl;
      result.scale(&celliptic, 1);
    }
  }
  
#endif





  // release allocated memory
  if (f) delete f;
  if (ft) delete ft;
  
  return 0;
}

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
	       const int jmax,
	       bool verbose=false)
{
  typedef typename ELLIPTIC::Index Index;
  typedef InfiniteVector<double,Index> V; // without "typename"...

  const unsigned int stages = row_stage_equation->row_method_->A.row_dimension(); // for readability

  std::list<V> Dalpha_uis; // will hold D_alpha u_i, i=1,...,s

  // solve the stage equations sequentially

  const double alpha = 1./(h*row_stage_equation->row_method_->C(0,0)); // we assume that the gamma_{i,i} are all the same
  if (verbose) cout << "increment: alpha=" << alpha << endl;
  row_stage_equation->set_alpha(alpha);

  for (unsigned int i(0); i < stages; i++) {  
    // setup right-hand side of i-th stage equation
    if (verbose) cout << "increment: setup rhs of stage #" << i << ":" << endl;
    row_stage_equation->setup_rhs(i, tolerance/(4*stages), t_n, h, D_un, Dalpha_uis, jmax);
    if (verbose) cout << "  ... done" << endl;

    // solve i-th stage equation
    V result;
#if 0
    if (verbose) cout << "increment: call CDD1_SOLVE()..." << endl;
    CDD1_SOLVE(*row_stage_equation, tolerance/(4*stages),
	       D_un, // guess
	       result,
 	       1.0,
 	       1.0,
	       jmax); // B_alpha*D_alpha x = D_alpha^{-1}y
    if (verbose) cout << "increment: CDD1_SOLVE() done" << endl;
#else
    GALERKIN_SOLVE(*row_stage_equation,
		   D_un,
		   result,
		   jmax);
#endif
    
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

/*!
  given an elliptic operator, u_0 and a stage equation object,
  solve u'(t)=Au(t)+f(t) adaptively
*/
template <class ELLIPTIC>
void solve_IVP(const ELLIPTIC* elliptic,
	       const InfiniteVector<double,typename ELLIPTIC::Index>& u0,
	       ROWStageEquation<ELLIPTIC>* row_stage_equation,
 	       const double T,
 	       const double atol,
	       const double rtol,
 	       const double q,
 	       const double tau_max,
 	       IVPSolution<InfiniteVector<double,typename ELLIPTIC::Index> >& result,
	       const int jmax)
{
  typedef typename ELLIPTIC::Index Index;
  typedef InfiniteVector<double,Index> V; // without "typename"...

  // IVP solver a la Hairer/Wanner
  result.t.clear();
  result.u.clear();
  
  double t_m = 0;
  V u_m(u0);
  
  unsigned int m = 0;
  result.t.push_back(t_m);
  result.u.push_back(u_m);
  
  cout << "t_{" << m << "}=" << t_m << " accepted!" << endl;
  
  m++;
    
  const double rho = 0.8; // overall safety factor

  // guess the initial time stepsize (cf. Hairer/Wanner, p. 169)
//   double u0_norm = linfty_norm(u0);
//   VECTOR yp0(u0);
//     ivp->evaluate_f(0, ivp->u0, atol, yp0); // initial slope
//     double ft0u0_norm = linfty_norm(yp0);
//     double tau0 = (u0_norm < 1e-5 || ft0u0_norm < 1e-5)
//       ? 1e-6 : 1e-2*u0_norm/ft0u0_norm;

  double tau0 = 1e-6; // pessimistic

//     // TODO: approximate ypp
    
  double tau_m = std::min(tau_max, 100*tau0);
    
  V u_mplus1, error_estimate;
  bool done = false;

  // parameters for step size selection;
  const double fac1 = 5.0;
  const double fac2 = 1./q;
  double fac = 0, facgus = 0, hacc = 0, erracc = 0; // will be set in the loop

  while (!done)
    {
      // jump to T if within 10% of T-t_m
      if (1.1*tau_m >= fabs(T-t_m) || t_m+tau_m >= T)
	{
	  tau_m = T-t_m;
	  done = true; // we would be "done" if the step is accepted
	}

      // try to advance one step
      u_m.scale(elliptic, 1); // u_m <- Du_m, coeffs w.r.t. D^{-1}Psi
      increment(elliptic, row_stage_equation,
		1e-2, t_m, tau_m, u_m, u_mplus1, error_estimate, jmax);
      u_m.scale(elliptic, -1); // u_m restored to L_2 coeffs
      u_mplus1.scale(elliptic, -1); // u_mplus1 now L_2 coeffs
//       error_estimate.scale(elliptic, -1); // error_estimate in L_2
      
      // compute an error estimator
// 	const double u_m_norm = linfty_norm(u_m);
// 	const double u_mplus1_norm = linfty_norm(u_mplus1);
// 	const double maxnorm = std::max(u_m_norm, u_mplus1_norm);
//  	double errest = 0;
//  	for (typename VECTOR::const_iterator it(error_estimate.begin());
// 	     it != error_estimate.end(); ++it)
// 	  errest += ((*it * *it)/((atol+maxnorm*rtol)*(atol+maxnorm*rtol)));
	
// 	errest = sqrt(errest / error_estimate.size());
//  	errest = linfty_norm(error_estimate) / (atol+maxnorm*rtol); // not good

      double errest = error_estimate.wrmsqr_norm(atol, rtol, u_m, u_mplus1);

      // estimate new stepsize
      fac = std::max(fac2, std::min(fac1, pow(errest, 1./row_stage_equation->row_method_->order()) / rho));
      double tau_new = tau_m / fac;
	
      if (errest <= 1)
	{
	  // accept the time step

	  cout << "t_{" << m << "}=" << t_m+tau_m << " accepted!" << endl;

	  t_m += tau_m;
	  result.t.push_back(t_m);
	  u_m = u_mplus1;
	  result.u.push_back(u_m);

	  // predictive controller of Gustafsson
	  if (m >= 2) {
	    facgus = (hacc/tau_m) * pow(errest*errest/erracc, 1./row_stage_equation->row_method_->order()) / rho;
	    facgus = std::max(fac2, std::min(fac1, facgus));
	    fac = std::max(fac, facgus);
	  }
	  hacc = tau_m;
	  erracc = std::max(1e-2, errest);

	  tau_m = tau_new;
	    
	  m++;
	}
      else
	{
	  // reject the time step

	  cout << "t_{" << m << "}=" << t_m+tau_m << " rejected!"
	       << " (errest=" << errest << ")" << endl;
	    
	  tau_m = tau_new;
	    
	  done = false;
	}
    }
}

/*!
  nonadaptive Galerkin solver, up to a fixed maximal level
*/
template <class PROBLEM>
void GALERKIN_SOLVE(const PROBLEM& P,
		    const InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& guess,
		    InfiniteVector<double, typename PROBLEM::WaveletBasis::Index>& u_epsilon,
		    const int jmax)
{
  typedef typename PROBLEM::WaveletBasis::Index Index;

  // setup index set Lambda
  set<Index> Lambda;
  for (Index lambda = P.basis().first_generator(P.basis().j0());; ++lambda) {
    Lambda.insert(lambda);
    if (lambda == P.basis().last_wavelet(P.basis().j0())) break;
  }

  // setup stiffess matrix
  SparseMatrix<double> A_Lambda;
  setup_stiffness_matrix(P, Lambda, A_Lambda);

  // setup righthand side
  Vector<double> F_Lambda;
  setup_righthand_side(P, Lambda, F_Lambda);
  
  // solve Galerkin system with CG
  Vector<double> xk(Lambda.size());
  unsigned int iterations = 0;
  const double eta = 1e-5;
  CG(A_Lambda, F_Lambda, xk, eta, 150, iterations);
  
  unsigned int id = 0;
  for (typename set<Index>::const_iterator it = Lambda.begin(), itend = Lambda.end();
       it != itend; ++it, ++id)
    u_epsilon.set_coefficient(*it, xk[id]);
}

int main()
{
  cout << "Testing ROWStageEquation..." << endl;

  // setup elliptic operator -Delta
  Poisson1D poisson_bvp;

  // setup 1D wavelet basis
  const int d  = 3;
  const int dT = 3;
//   typedef DSBasis<d,dT> Basis; string basis_str = "DS";
  typedef PBasis<d,dT> Basis; string basis_str = "P";

  const int jmax = 9;

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
  ROWMethod<Vector<double> > row_method(WMethod<Vector<double> >::ROS2); string scheme_str="ROS2";
//   ROWMethod<Vector<double> > row_method(WMethod<Vector<double> >::ROS3); string scheme_str="ROS3";
//   ROWMethod<Vector<double> > row_method(WMethod<Vector<double> >::RODASP); string scheme_str="RODASP";


//     ROWMethod<V> row_adaptive(WMethod<V>::ROS2);
//     ROWMethod<V> row_adaptive(WMethod<V>::ROS3P);
//     ROWMethod<V> row_adaptive(WMethod<V>::ROS3Pw);
//     ROWMethod<V> row_adaptive(WMethod<V>::ROSI2P2);
//     ROWMethod<V> row_adaptive(WMethod<V>::ROS3);
//     ROWMethod<V> row_adaptive(WMethod<V>::GRK4T);
//     ROWMethod<V> row_adaptive(WMethod<V>::ROWDA3);
//     ROWMethod<V> row_adaptive(WMethod<V>::RODASP);



  //
  //
  // setup the ROW stage equation object
  ROWStageEquation<CachedProblem<EllipticEquation> > row_stage_equation(&row_method, &celliptic, f, ft);


#if 0
  // nonadaptive solution with constant temporal stepsize h

  cout << "* testing approximation with constant stepsizes..." << endl;

  V temp, error_estimate, result;
  IVPSolution<V> results;

  const int expo1 = 9;
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
      temp = result;
      results.t.push_back(i*h);
      result.scale(&celliptic, -1); // switch to the L_2 coeffs
      results.u.push_back(result);

      ostringstream output_filename;
      output_filename << "u" << i << ".m";
      resultstream.open(output_filename.str().c_str());
      SampledMapping<1> ui_plot(evaluate(celliptic.basis(), result, true, resolution));
      ui_plot.matlab_output(resultstream);
      resultstream.close();

      cout << "---------------- after increment() -----------------------" << endl;

      V uexact_coeffs;
      uexact.set_time(i*h);
      expand(&uexact, celliptic.basis(), false, jmax, uexact_coeffs);
      cout << "  ell_2 error at t=" << i*h << ": " << l2_norm(result - uexact_coeffs) << endl;
    }
  }
  
#else

  const double T = 1.0;
  const double q = 10.0;
  const double tau_max = 1.0;
  
  std::list<double> numberofsteps;
  std::list<double> errors;
  std::list<double> wallclocktimes;

  cout << "* testing scheme " << scheme_str << " (adaptive, several tolerances)..." << endl;
  for (int expo = 6; expo <= 18; expo++) { // 2^{-6}=0.015625, 2^{-8}=3.9e-3, 2^{-10}=9.77e-4
    const double TOL = ldexp(1.0, -expo);
//     const double TOL = pow(10.0, -(double)expo);
    
    IVPSolution<V> result_adaptive;
    
    cout << "  TOL=" << TOL << endl;

    clock_t tstart =  clock();
    solve_IVP(&celliptic, u0, &row_stage_equation, T,
	      TOL, 0, q, tau_max, result_adaptive, jmax);
    clock_t tend = clock();

#if _TESTCASE == 1
    // compute maximal ell_2 error of the coefficients
    double errhelp = 0;
    std::list<double>::const_iterator ti(result_adaptive.t.begin());
    for (std::list<V>::const_iterator ui(result_adaptive.u.begin());
	 ui != result_adaptive.u.end(); ++ui, ++ti) {
      V uexact_coeffs;
      uexact.set_time(*ti);
      expand(&uexact, celliptic.basis(), false, jmax, uexact_coeffs);
      cout << "  ell_2 error at t=" << *ti << ": " << l2_norm(*ui - uexact_coeffs) << endl;
      errhelp = std::max(errhelp, l2_norm(*ui - uexact_coeffs));
    }
    errors.push_back(errhelp);
    numberofsteps.push_back(result_adaptive.t.size());
    wallclocktimes.push_back((double)(tend-tstart)/ (double)CLOCKS_PER_SEC);
#endif

    ostringstream filename;
    filename << basis_str << "_" << scheme_str
	     << "_case" << _TESTCASE
	     << "_d" << d <<  "_dt" << dT
	     << "_jmax" << jmax
	     << ".m";
    std::ofstream resultstream(filename.str().c_str());
    
    resultstream << "N_errors=[";
    for (std::list<double>::const_iterator it = errors.begin();
	 it != errors.end(); ++it) {
      resultstream << log10(*it);
      if (it != errors.end())
	resultstream << " ";
    }
    resultstream << "];" << endl;
    
    resultstream << "N=[";
    for (std::list<double>::const_iterator it = numberofsteps.begin();
	 it != numberofsteps.end(); ++it) {
      resultstream << log10(*it);
      if (it != numberofsteps.end())
	resultstream << " ";
    }
    resultstream << "];" << endl;
        
    resultstream << "time_errors=[";
    for (std::list<double>::const_iterator it = errors.begin();
	 it != errors.end(); ++it) {
      resultstream << log10(*it);
      if (it != errors.end())
	resultstream << " ";
    }
    resultstream << "];" << endl;
    
    resultstream << "times=[";
    for (std::list<double>::const_iterator it = wallclocktimes.begin();
	 it != wallclocktimes.end(); ++it) {
      resultstream << log10(*it);
      if (it != wallclocktimes.end())
	resultstream << " ";
    }
    resultstream << "];" << endl;
    
  }
  
#endif
  
  // release allocated memory
  if (f) delete f;
  if (ft) delete ft;
  
  return 0;
}

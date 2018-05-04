//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!  The Kacanov-type iteration method for the p-Poisson problem for 1<p<2 from [Diening, Fornasier, Wank: 'A Relaxed Kacanov Iteration for the p-Poisson Problem' (2017)]
//!
//!  Tests on the cube (the adaptive wavelet-Galerkin method (CDD1_SOLVE) is used for the solution of the arising linear subproblems)
//!
//!  Note: For the storage of the computed data (approximation errors, plots, parameter values, etc.), the subfolder structure './results/plots' must exist in the folder of the executable
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#define CP_USE_PRECOMPUTED_NORM 1           // Cached problem: use precomputed eigenvalue bounds of stiffness matrix
#define CEPP_USE_ALTERNATIVE_F 1           // cube_equation_pPoisson, member f: use alternative representation of f
#define _DIM 2
#define P_POISSON

#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <numerics/bvp.h>
#include <numerics/iteratsolv.h>

#include <interval/i_index.h>
#include <interval/p_basis.h>

#include <cube/cube_indexplot.h>

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 1

#include <cube/cube_basis.h>
#include <galerkin/cached_problem.h>
#include <galerkin/cube_equation_pPoisson.h>

#define _WAVELETTL_CDD1_VERBOSITY 1
#include <adaptive/cdd1.h>
#include <galerkin/p_poisson_data_cube.h>

using namespace std;
using namespace WaveletTL;
using namespace MathTL;

//! external variables for time measurement
double time_consumption_of_a;
double time_consumption_of_NGROW;
double time_consumption_of_GALERKIN;
double time_consumption_of_NRESIDUAL;
double time_consumption_of_coeff_func;
double time_consumption_of_calc_deriv;
double time_consumption_of_calc_deriv_outer;
double time_consumption_of_compute_rhs;
double time_consumption_of_eval_u_eps;
double time_consumption_of_eval_deriv_u_eps;

//! external variables for information about cache usage (cache for entries of stiffness matrix)
int number_of_entries_computed;
int number_of_entries_from_cache;



//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!   classes ('Exact_Solution', 'Deriv_Exact_Solution', 'TestRHS', 'TestRHS_AR') for some test problems for the p-Poisson equation on the cube with homogeneous Dirichlet b.c.'s:
//!
//!   1: A smooth, polynomial solution: u(x,y) = x*(x-1)*y*(y-1)
//!   2: The 'peak': u(x,y) = max{0.5 - |(x,y) - (0.5, 0.5)|, 0}
//!
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! the exact solution $u$ of the p-Poisson equation
template <unsigned int N>
class Exact_Solution
  : public Function<2,double>
{
public:
  virtual ~Exact_Solution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const
  {
    switch(N)
    {
    case 1:
      return p[0]*(p[0]-1)*p[1]*(p[1]-1);
      break;
    case 2:
      return max((0.5 - sqrt( (p[0] - 0.5) * (p[0] - 0.5) + (p[1] - 0.5) * (p[1] - 0.5) )), 0.0);
      break;
    case 3:
      return 0.0;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const {
    values[0] = value(p);
  }
};

//! derivative of the exact solution $u$
template <unsigned int N>
class Deriv_Exact_Solution
  : public Function<2,double>
{
public:
  virtual ~Deriv_Exact_Solution() {};
  double value(const Point<2>& p, const unsigned int component = 0) const
  {
    switch(N)
    {
    case 1:
      return sqrt( (2*p[0] - 1) * (p[1]*p[1] - p[1]) * (2*p[0] - 1) * (p[1]*p[1] - p[1]) + (2*p[1] - 1) * (p[0]*p[0] - p[0]) * (2*p[1] - 1) * (p[0]*p[0] - p[0]) );
      break;
    case 2:
      {
        double r = sqrt( (p[0] - 0.5) * (p[0] - 0.5) + (p[1] - 0.5) * (p[1] - 0.5) );
        if ( r < 0.5 )
        {
          return 1.0;
        }
        else
        {
          return 0.0;
        }
      }
      break;
    case 3:
      return 0.0;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const
  {
    switch(N)
    {
    case 1:
      {
        values[0] = (2*p[0] - 1) * (p[1]*p[1] - p[1]);
        values[1] = (2*p[1] - 1) * (p[0]*p[0] - p[0]);
      }
      break;
    case 2:
      {
        double r = sqrt( (p[0] - 0.5) * (p[0] - 0.5) + (p[1] - 0.5) * (p[1] - 0.5) );
        if (r < 0.5)
        {
          if (r > 1.0e-14)
          {
            values[0] = (0.5 - p[0]) / r;
            values[1] = (0.5 - p[1]) / r;
          }
          else
          {
            values[0] = 0.0;
            values[1] = 0.0;
          }
        }
        else
        {
          values[0] = 0.0;
          values[1] = 0.0;
        }
      }
      break;
    case 3:
      {
      values[0] = 0.0;
      values[1] = 0.0;
      }
      break;
    }
  }
};

/*!
  Right-hand side (rhs) $f$ of the p-Poisson equation:
  f := -Delta_p u(x,y) = -div(|nabla_u|^(p-2) * nabla_u)
*/
template <unsigned int N>
class TestRHS
  : public Function<2,double>
{
public:
  virtual ~TestRHS() {};
  double param_p;

  double value(const Point<2>& p, const unsigned int component = 0) const
  {
    switch(N)
    {
      case 1:
      {
        double u_x = (2*p[0]-1)*p[1]*(p[1]-1);
        double u_y = (2*p[1]-1)*p[0]*(p[0]-1);
        double u_xx = 2*p[1]*(p[1]-1);
        double u_yy = 2*p[0]*(p[0]-1);
        double u_xy = (2*p[0]-1)*(2*p[1]-1);
        double nabla_u_square = u_x*u_x + u_y*u_y;
        double r;

        r = (param_p-2) * pow(nabla_u_square, (param_p-4)*0.5) * (u_x*u_x*u_xx + u_y*u_y*u_yy + 2*u_x*u_y*u_xy)
          + pow(nabla_u_square, (param_p-2)*0.5) * (u_xx + u_yy);
        return -r;
      }
      break;
      case 2:
      {
        double r = sqrt( (p[0] - 0.5) * (p[0] - 0.5) + (p[1] - 0.5) * (p[1] - 0.5) );
        if (r < 0.5)
        {
          if (r > 1.0e-14)
          {
            return (1.0 / r);
          }
          else
          {
            return 1.0e14;
          }
        }
        else
        {
          return 0.0;
        }
      }
      break;
      case 3:
        return 0.0;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const
  {
    values[0] = value(p);
  }
};

/*!
  Alternative representation of the rhs $f$:
  returns $\abs{\nabla u}^{p-2} \nabla u$, which can be used for, e.g., the computation of $f(v) = \int_{\Omega} \abs{\nabla u}^{p-2} \nabla u * \nabla v dx$
*/
template <unsigned int N>
class TestRHS_AR
  : public Function<2,double>
{
public:
  virtual ~TestRHS_AR() {};
  double param_p;

  double value(const Point<2>& p, const unsigned int component = 0) const
  {
    switch(N)
    {
      case 1:
      {
        return 0.0;
      }
      break;
      case 2:
      {
        return 0.0;
      }
      break;
      case 3:
        return 0.0;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const
  {
    switch(N)
    {
      case 1:
      {
        values[0] = 0.0;
        values[1] = 0.0;
      }
      break;
      case 2:
      {
        double r = sqrt( (p[0] - 0.5) * (p[0] - 0.5) + (p[1] - 0.5) * (p[1] - 0.5) );
        if (r < 0.5)
        {
          if (r > 1.0e-14)
          {
            values[0] = (0.5 - p[0]) / r;
            values[1] = (0.5 - p[1]) / r;
          }
          else
          {
            values[0] = 0.0;
            values[1] = 0.0;
          }
        }
        else
        {
          values[0] = 0.0;
          values[1] = 0.0;
        }
      }
      break;
      case 3:
      {
        values[0] = 0.0;
        values[1] = 0.0;
      }
      break;

    }
  }
};

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!  class 'Coeff_pPoisson' models the coefficient function of an elliptic PDE with nonconstant coefficients:
//!  tailored version for the (relaxed) p-Poisson operator and representation of a(x) in frame coordinates: member 'value' returns
//!  ( min{ max{cutoff_small, |\nabla a(x)|}, cutoff_big} )^{p-2},   where 0 < cutoff_small < cutoff_big < infty are the relaxation parameters
//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
template <class IBASIS, unsigned int DIM>
class Coeff_pPoisson
  : public Function<DIM,double>
{
public:
  virtual ~Coeff_pPoisson() {};

  //! Constructor from bc's
  Coeff_pPoisson( const FixedArray1D<bool,2*DIM>& bc )
    : _basis(bc)
  {
  }

  //! Constructor from bc's
  Coeff_pPoisson( const FixedArray1D<int,2*DIM>& bc )
    : _basis(bc)
  {
  }

  double param_p;

  //! point evaluation
  double value(const Point<DIM>& p, const unsigned int component = 0) const
  {

    FixedArray1D<double,DIM> deriv_values;
    double nabla_square_norm(0.);
    double r;

//    evaluate(_basis, a_coeffs, p, deriv_values);  //! compute derivative of a_coeffs in point p and store in deriv_values
//    evaluate_tuned(_basis, a_coeffs, p, deriv_values);  //! tuned Version
//    evaluate_tuned_dim2(_basis, a_coeffs, p, deriv_values);  //! tuned Version (only for dimension=2)

////! for sanity-check for evaluate
//FixedArray1D<double,DIM> deriv_values2;
//double diff0, diff1;
////! for sanity-check for evaluate ENDE

    // measure time for function call
    clock_t begin_EVALUATE_a = clock();

    evaluate_tuned_dim2_lin_primbs(_basis, a_coeffs, p, deriv_values);  //! tuned Version (only for dimension=2, linear (2,2) Primbs basis with homog. bc's)

    // measure time for function call "evaluate_tuned_dim2_2"
    clock_t end_EVALUATE_a = clock();
    double elapsed_secs = double(end_EVALUATE_a - begin_EVALUATE_a) / CLOCKS_PER_SEC;
    time_consumption_of_coeff_func += elapsed_secs;

////! sanity-check for evaluate
//    diff0 = deriv_values[0] - deriv_values2[0];
//    diff1 = deriv_values[1] - deriv_values2[1];
//
//    cout << setprecision(20);
//    if (abs(diff0) > 1e-15)
//    {
//      cout << "ERROR: deriv_values mismatch in x!!" << endl;
//      cout << deriv_values[0] << " and " << deriv_values2[0] << endl;
//      exit(1);
//    }
//
//    if (abs(diff1) > 1e-15)
//    {
//      cout << "ERROR: deriv_values mismatch in y!!" << endl;
//      cout << deriv_values[1] << " and " << deriv_values2[1] << endl;
//      exit(1);
//    }
//    cout << setprecision(6);
////! sanity-check ENDE

    for (unsigned int i = 0; i < DIM; i++)  //! compute l_2-norm of derivative
    {
      nabla_square_norm += (deriv_values[i]*deriv_values[i]);
    }

    r = sqrt(nabla_square_norm);

    if (r < cutoff_small)
      r = cutoff_small;
    else if (r > cutoff_big)
      r = cutoff_big;

    return pow(r, param_p - 2.0);

  }

  void vector_value(const Point<DIM>& p, Vector<double>& values) const
  {
    values[0] = value(p);
  }

  void setcoeff(const double epsilon_lb, const double epsilon_ub, const InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index> coeffs)
  {
    cutoff_small = epsilon_lb;
    cutoff_big = epsilon_ub;
    a_coeffs.clear();
    a_coeffs += coeffs;
  }

  //! data
  double cutoff_small, cutoff_big;
  //FixedArray1D<double,DIM> deriv_values2;
  InfiniteVector<double, typename CubeBasis<IBASIS,DIM>::Index> a_coeffs;
  CubeBasis<IBASIS,DIM> _basis;

};

//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


int main()
{
  cout << "Kacanov-type iteration method for the p-Poisson problem on the cube (the adaptive wavelet-Galerkin method (CDD1_SOLVE) is used for the solution of the linear subproblems) ..." << endl;

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                   Parameters and flags
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const int PARAM_EXAMPLE = 2;                             // Encodes the p-Poisson equation to be solved (all on the cube with homog. Dirichlet bc's). The corresponding right-hand
                                                           // side is defined in class 'TestRHS' (and TestRHS_AR resp.), the explicit solution (if known) in class 'Exact_Solution'.

  double param_p = 1.5;                      // the parameter 'p' of the p-Poisson problem

  const int d  = 2;
  const int dT = 2;

  //const int MAX_LEV_WAV_ALL = 7;                           // the maximal level of the wavelets
const int MAX_LEV_WAV_ALL = 4;  // few levels for testing (speedup)

  const int max_lev_wav_basis = MAX_LEV_WAV_ALL;           // the maximal level of the wavelet basis (CUBEBASIS)
  const int max_lev_wav_rhs = MAX_LEV_WAV_ALL;             // the maximal level of the wavelets considered by 'compute_rhs()' for precomputation of the rhs
  const int max_lev_wav_CDD1 = MAX_LEV_WAV_ALL;            // the maximal level of the wavelets considered by 'CDD1_SOLVE'
  const int max_lev_wav_u_exact = MAX_LEV_WAV_ALL;         // the maximal level of the wavelets considered for the expansion of the exact solution u_exact

  //const int max_kacanov_iterations = 15;     // maximal number of Kacanov iterations. in each iteration we solve a lin. elliptic PDE with nonconst. coeff.
const int max_kacanov_iterations = 2; // few iterations for testing (speedup)

  const int min_res_quad_a = 0;              // parameter for the evaluation of the bilinear form a(,). maximal side length (=2^-min_resolution_quadrature) of the cubes
                                             // used for the composite Gauss quadrature
  const int doe_quad_a = 3;                  // degree of exactness of the Gauss quadrature rule used in 'a' (Note: only odd doe possible. For even value, doe will be value-1 !!)

  //const int min_res_quad_f = 8;              // parameter for the evaluation of the functional f(). maximal side length (=2^-min_resolution_quadrature_f) of the cubes
                                               // used for the composite Gauss quadrature. Set value > 0 if rhs is not very smooth.
const int min_res_quad_f = 3; // bad resolution for testing (speedup)

  const int doe_quad_f = 3;                  // degree of exactness of the Gauss quadrature rule used in 'f' (Note: only odd doe possible. For even value, doe will be value-1 !!)

  double epsilon_stabilization = 0.001;         //  the stabilization (relaxation) parameter for the Kacanov iteration

  //double epsilon_CDD1 = 0.001;                // The target $\ell_2$-accuracy for the linear subproblems. 'CDD1_SOLVE' terminates if ||u - u_epsilon||_{l_2} < epsilon_CDD1
double epsilon_CDD1 = 0.1; // bad accuracy for testing (speedup)

  double gamma_CDD1 = 0.85;

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                     typedef's
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  // the wavelet basis (on the square [0,1]^2) that we use
  typedef PBasis<d,dT> Basis1D;             // Primbs basis on the intervall [0,1]
  typedef CubeBasis<Basis1D,2> Basis;       // tensor product wavelet basis on [0,1]^2
  typedef Basis::Index Index;

  typedef CubeEquationpPoisson<Basis1D,2,Basis> Problem_pPoisson;

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                 Variable declarations
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  //! the new and old approximant
  InfiniteVector<double, Index> u_epsilon;
  InfiniteVector<double, Index> u_epsilon_old;

  //! bc's
  FixedArray1D<bool,4> bc;
  bc[0] = bc[1] = bc[2] = bc[3] = true;     // i.e., homogeneous Dirichlet b.c.'s

  //! right-hand side f: \Omega -> \R (as pointwise function), i.e.: rhs_{\lambda} = \int_{\Omega} f(x) * Psi_{\lambda}(x) dx
  TestRHS<PARAM_EXAMPLE> rhs;
  rhs.param_p = param_p;

  //! right-hand side f_ar: \Omega -> \R^d (alternative representation of the rhs functional), i.e.: rhs_{\lambda} = \int_{\Omega} f_ar(x) * \nabla Psi_{\lambda}(x) dx
  TestRHS_AR<PARAM_EXAMPLE> rhs_ar;
  rhs_ar.param_p = param_p;

  //! exact solution
  Exact_Solution<PARAM_EXAMPLE> u_exact;                            // exact solution as a function
  Deriv_Exact_Solution<PARAM_EXAMPLE> deriv_u_exact;                // derivative of the exact solution as a function

  //! coefficient a
  Coeff_pPoisson<Basis1D,2> mein_coeff_pP(bc);
  mein_coeff_pP.param_p = param_p;
  mein_coeff_pP.setcoeff(epsilon_stabilization, 1.0/epsilon_stabilization, InfiniteVector<double, Index>());

  //! setup bvp
  PoissonBVP_Coeff<2> poisson_coeff(&mein_coeff_pP, &rhs, &rhs_ar);

  //! setup problem
  clock_t begin_compute_rhs = clock();

  Problem_pPoisson problem(&poisson_coeff, bc, param_p, max_lev_wav_basis, max_lev_wav_rhs, 
          doe_quad_a, min_res_quad_a, doe_quad_f, min_res_quad_f);

  clock_t end_compute_rhs = clock();
  time_consumption_of_compute_rhs = double(end_compute_rhs - begin_compute_rhs) / CLOCKS_PER_SEC;



  //! setup cached problem
#if CP_USE_PRECOMPUTED_NORM
 //  CachedProblem<Problem_pPoisson> cproblem(&problem, 2.6, 12.0);  // initialization with precomputed PBasis eigenvalue bounds (d=dt=2) (EXAMPLE 1)
   CachedProblem<Problem_pPoisson> cproblem(&problem, 2.9, 12.0);  // initialization with precomputed PBasis eigenvalue bounds (d=dt=2) (EXAMPLE 2)
#else
   CachedProblem<Problem_pPoisson> cproblem(&problem);
#endif // CP_USE_PRECOMPUTED_NORM


  // instantiate class 'pPoissonData' for computation and storage of: the approximation error (in various norms), plots, parameters, energy...
  pPoissonData<CachedProblem<Problem_pPoisson> > p_poisson_data(&cproblem, &u_exact, &deriv_u_exact);

  /* compute:   - the wavelet expansion of the exact solution (up to wavelet level 'max_lev_wav_u_exact'),
                - the exact solution as a sampled mapping and as a matrix,
                - a Matlab plot of the exact solution,
                - the energy of the exact solution
  */
  p_poisson_data.precompute_exact_solution(max_lev_wav_u_exact, 3, 8);

//! +++++++++++++++++++++++++++++++++++++++++++++++++++++     *END* of Variable declarations     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                  Start of the Kacanov iteration for solving p-Poisson
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  for (int i = 0; i < max_kacanov_iterations; i++)
  {

    time_consumption_of_a = 0.;
    time_consumption_of_NGROW = 0.;
    time_consumption_of_GALERKIN = 0.;
    time_consumption_of_NRESIDUAL = 0.;
    time_consumption_of_coeff_func = 0.;
    time_consumption_of_calc_deriv = 0.;
    time_consumption_of_calc_deriv_outer = 0.;
    time_consumption_of_eval_u_eps = 0.;
    time_consumption_of_eval_deriv_u_eps = 0.;


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! compute next iterate u_epsilon

    clock_t begin_CDD1 = clock();

    CDD1_SOLVE(cproblem, epsilon_CDD1, u_epsilon_old, u_epsilon, max_lev_wav_CDD1, CDD1, gamma_CDD1);


    //! scale u_epsilon
    u_epsilon.scale(&cproblem, -1);


    clock_t end_CDD1 = clock();
    double time_consumption_of_CDD1 = double(end_CDD1 - begin_CDD1) / CLOCKS_PER_SEC;


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! save computed data

    clock_t begin_compute_data = clock();

    cout << "'p_poisson_data': compute_data .." << endl;

    p_poisson_data.compute_data(u_epsilon, param_p, 3, 8, true, epsilon_stabilization, 1.0/epsilon_stabilization);

    cout << "'p_poisson_data': .. finished!" << endl;

    clock_t end_compute_data = clock();
    double time_consumption_of_compute_data = double(end_compute_data - begin_compute_data) / CLOCKS_PER_SEC;

    p_poisson_data.Data[p_poisson_data.iterations() - 1].epsilon_CDD1 = epsilon_CDD1;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].epsilon_stabilization = epsilon_stabilization;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].max_lev_wav_CDD1 = max_lev_wav_CDD1;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].max_lev_wav_rhs = max_lev_wav_rhs;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].max_lev_wav_u_exact = max_lev_wav_u_exact;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].min_res_quad_a = min_res_quad_a;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].doe_quad_a = doe_quad_a;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].min_res_quad_f = min_res_quad_f;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].doe_quad_f = doe_quad_f;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].gamma_CDD1 = gamma_CDD1;

    cout << "iteration " << i << ", max. absolute error of sampled approximation = " << p_poisson_data.Data[i].abs_error_sm_max_nrm << endl;
    cout << "iteration " << i << ", max. relative error of sampled approximation = " << p_poisson_data.Data[i].rel_error_sm_max_nrm << endl;
    cout << "iteration " << i << ", \ell_2 error of u_epsilon = " << p_poisson_data.Data[i].error_l2 << endl;
    cout << "iteration " << i << ", \ell_infty error of u_epsilon = " << p_poisson_data.Data[i].error_l_infty << endl;


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! output time measurement

    cout << endl;
    cout << "time for compute_rhs, i = " << i << " : " << time_consumption_of_compute_rhs << endl;
    cout << "time for CDD1, i = " << i << " : " << time_consumption_of_CDD1 << endl;
    cout << "time for NGROW, i = " << i << " : " << time_consumption_of_NGROW << endl;
    cout << "time for NRESIDUAL, i = " << i << " : " << time_consumption_of_NRESIDUAL << endl;
    cout << "time for GALERKIN, i = " << i << " : " << time_consumption_of_GALERKIN << endl;
    cout << "time for a(), i = " << i << " : " << time_consumption_of_a << endl;
    cout << "time for coeff_func, i = " << i << " : " << time_consumption_of_coeff_func << endl;
    cout << "time for calc_deriv + outer loop, i = " << i << " : " << time_consumption_of_calc_deriv_outer << endl;
    cout << "time for calc_deriv, i = " << i << " : " << time_consumption_of_calc_deriv << endl;

    cout << "time for compute_data, i = " << i << " : " << time_consumption_of_compute_data << endl;
    cout << "time for eval_u_eps, i = " << i << " : " << time_consumption_of_eval_u_eps << endl;
    cout << "time for eval_deriv_u_eps, i = " << i << " : " << time_consumption_of_eval_deriv_u_eps << endl;

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Perform several updates before next iteration


    if (i == (max_kacanov_iterations-1))
    {
      break;
    }

    // set minimal resolution of composite gauss quadrature to: (highest wavelet level of u_epsilon) + 1
    std::set<Index> supp_u;
    u_epsilon.support(supp_u);
    std::set<Index>::const_iterator it = supp_u.end();
    Index dummy_index;

    if (supp_u.empty())
    {
      //j = frame->j0() + 1;
      problem.min_res_quadrature_a = 3 + 1;
    }

    else
    {
      --it;
      dummy_index = *it;
      problem.min_res_quadrature_a = dummy_index.j() + 1;
    }


    // some output..
    if (supp_u.empty())
    {
      cout << "p_Poisson: last = null" << " level = " << 3 << endl;
    }
    else
    {
      cout << "p_Poisson: last = " << dummy_index << " level = " << dummy_index.j() << endl;
    }
    cout << "p_Poisson: i=" << i << ", min_res_quadrature_a: " << problem.min_res_quadrature_a << endl;

    // update stabilization parameter
    //epsilon_stabilization *= 0.5;

    // update coefficient a
    mein_coeff_pP.setcoeff(epsilon_stabilization, 1.0/epsilon_stabilization, u_epsilon);    // update coefficient function
    poisson_coeff.set_a(&mein_coeff_pP);             // upadate bvp

#if !CP_USE_PRECOMPUTED_NORM
    problem.set_normA(0.);
    problem.set_normAinv(0.);
    cproblem.set_normA(0.);
    cproblem.set_normAinv(0.);
#endif

    cproblem.clear_cache();                    // clear the entries cache for A
    problem.clear_coeff_cache();               // clear cache for the coefficient function

  clock_t begin_compute_rhs = clock();

    problem.set_bvp(&poisson_coeff);           // set bvp (-> invokes call of compute_rhs)

  clock_t end_compute_rhs = clock();
  time_consumption_of_compute_rhs = double(end_compute_rhs - begin_compute_rhs) / CLOCKS_PER_SEC;


    u_epsilon_old = u_epsilon;

  }

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!  *END* of the Kacanov iteration for solving p-Poisson                                   +
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//! some OUTPUT:

  cout << "Basis1D: primal_polynomial_degree: " << Basis1D::primal_polynomial_degree() << endl;

  //! save parameters, computed errors etc. in file. Output to terminal as well.
  std::ofstream ofs_pPoisson_data("./results/pPoisson_data.txt");                //! open file

  ofs_pPoisson_data << "p-Poisson equation: Results of example " << PARAM_EXAMPLE << " with p = " << param_p << endl << endl;
  for (int i = 0; i < p_poisson_data.iterations(); i++)
  {
    p_poisson_data.print_data(i, std::cout);
    p_poisson_data.print_data(i, ofs_pPoisson_data);
  }
  ofs_pPoisson_data.close();                                                     //! close file


  //! save some error plots */
  std::ofstream ofs_plot_error_energy("./results/error_energy.m");             //! open file
  std::ofstream ofs_plot_error_l2("./results/error_l2.m");                     //! open file
  std::ofstream ofs_plot_error_W1p("./results/error_W1p.m");                   //! open file
  std::ofstream ofs_plot_error_W1_infty("./results/error_W1_infty.m");         //! open file

  p_poisson_data.save_error_plots(ofs_plot_error_energy,                       //! save plots
                                  ofs_plot_error_l2,
                                  ofs_plot_error_W1p,
                                  ofs_plot_error_W1_infty);

  ofs_plot_error_energy.close();                                               //! close file
  ofs_plot_error_l2.close();                                                   //! close file
  ofs_plot_error_W1p.close();                                                  //! close file
  ofs_plot_error_W1_infty.close();                                             //! close file


//! rhs as matlab plot
//! rhs as sampled mapping
  Point<2, double> a(0.,0.), b(1., 1.);
  int number_gridpoints( 1 << 7);
  Grid<2> my_grid(a,b, number_gridpoints);
  SampledMapping<2> rhs_sm(my_grid, rhs);

  //! rhs as Matlab plot
  std::ofstream rhs_stream("./results/plots/rhs.m");
  rhs_sm.matlab_output(rhs_stream);
  rhs_stream.close();


  return 0;

}  //! *END* of main()



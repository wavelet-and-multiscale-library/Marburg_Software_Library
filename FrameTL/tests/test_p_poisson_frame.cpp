//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!  The Kacanov-type iteration method for the p-Poisson problem for 1<p<2 from [Diening, Fornasier, Wank: 'A Relaxed Kacanov Iteration for the p-Poisson Problem' (2017)]
//!
//!  Tests on the L-shaped domain ('Mult_Schw' is used for the solution of the arising linear subproblems)
//!
//!  Note: For the storage of the computed data (approximation errors, plots, parameter values, etc.), the subfolder structure './results/plots' must exist in the folder of the executable
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
#define CP_USE_PRECOMPUTED_NORM 1           // Cached problem: use precomputed eigenvalue bounds of stiffness matrix
#define PPE_USE_ALTERNATIVE_F 1             // class 'pPoissonEquation', member 'f': use alternative representation f_ar of f to compute rhs
#define PPDF_USE_ALTERNATIVE_F 1            // class 'pPoissonDataFrame', members 'precompute_exact_solution' and 'compute_data':
                                            // use alternative representation f_ar of f to compute second term of energy (energy2)

#define OVERLAP 1.

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 1
#define _WAVELETTL_CDD1_VERBOSITY 1

//#define PRECOMP_RHS
//#define PRECOMP_DIAG
#define COMPUTECONSTANTS 0

#define SPARSE
//#define FULL
#define TWO_D
#define _DIM 2
#define P_POISSON

#define EXAMPLE_SINGULARITY 1       // for PARAM_EXAMPLE = 2

#include <fstream>
#include <iostream>
#include <map>
#include <time.h>

#include <algebra/symmetric_matrix.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>

#include <numerics/bvp.h>
#include <numerics/iteratsolv.h>
#include <numerics/eigenvalues.h>
#include <numerics/corner_singularity.h>

#include <interval/i_index.h>
#include <interval/p_basis.h>

#include <cube/cube_indexplot.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>

#include <galerkin/cached_problem.h>
#include <galerkin/cube_equation_pPoisson.h>
#include <galerkin/galerkin_utils.h>

#include <adaptive/cdd1.h>
#include <pPoisson_frame_data.h>

#include <elliptic_equation.h>
#include <simple_elliptic_equation.h>
#include <pPoisson_equation.h>

#include <frame_evaluate.h>
#include <frame_support.h>
#include <frame_index.h>

#include <adaptive_multiplicative_Schwarz.h>
#include <error_H_scale.h>
//#include <multiplicative_Schwarz.h>


using std::cout;
using std::endl;

using FrameTL::FrameIndex;
using FrameTL::EllipticEquation;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
using FrameTL::pPoissonEquation;
using MathTL::EllipticBVP;
using MathTL::PoissonBVP;
using MathTL::ConstantFunction;
using MathTL::CornerSingularityRHS;
using MathTL::CornerSingularity;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;
using WaveletTL::CachedProblem;

using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


//! external variables for time measurement
double time_consumption_of_a;
double time_consumption_of_a_same_patches;
double time_consumption_of_a_different_patches;
double time_consumption_of_NGROW;
double time_consumption_of_GALERKIN;
double time_consumption_of_NRESIDUAL;
double time_consumption_of_coeff_func;
double time_consumption_of_calc_deriv;
double time_consumption_of_calc_deriv_outer;
double time_consumption_of_compute_rhs;
double time_consumption_of_compute_diagonal;
double time_consumption_of_eval_u_eps;
double time_consumption_of_eval_deriv_u_eps;

//! external variables for information about cache usage (cache for entries of stiffness matrix)
int number_of_entries_computed;
int number_of_entries_from_cache;



//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!   classes ('Exact_Solution', 'Deriv_Exact_Solution', 'TestRHS', 'TestRHS_AR') for some test problems for the p-Poisson equation on the L-domain with homogeneous Dirichlet b.c.'s:
//!
//!   1: A smooth, polynomial solution: $u(x,y) = x*(1-x^2)*y*(1-y^2)$
//!   2: The "singularity function": $u(r,theta) = zeta(r) * r^{2/3} * sin(2/3 * theta)$, where zeta is a smooth cut-off function and (r,theta) are polar coordinates w.r.t. the origin
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
      return (1 - (p[0]*p[0]) )*p[0]* (1 - (p[1]*p[1]) )*p[1];
      break;
    case 2:
    {
      double r = hypot(p[0],p[1]);

      if (r >= r1) return 0.0;

      double theta = atan2(p[1],p[0]);

      // shift theta to [0,2*pi]
      if (theta < 0) theta += 2.0 * M_PI;
      theta -= theta0 * M_PI;
      if (theta < 0) theta += 2.0 * M_PI;
      if (theta >= omega * M_PI) return 0.0;

      // compute zeta(r):
      double zeta = 0.0;
      if (r <= r0)
      {
        zeta = 1.0;
      }
      else
      {
        if (r >= r1)
        {
	      zeta = 0.0;
        }
        else
        {
          const double help1(r1-r);
          const double help0(r-r0);
	      zeta = exp(-1.0/(help1*help1))/(exp(-1.0/(help0*help0))+exp(-1.0/(help1*help1)));
        }
      }
      // END of compute zeta(r)

      return zeta * pow(r, 1.0/omega) * sin(theta/omega);
      break;
    }
    case 3:
      return 0.0;
      break;
    }
  }
  void vector_value(const Point<2>& p, Vector<double>& values) const
  {
    values[0] = value(p);
  }

  const static double r0 = 0.01;
  const static double r1 = 0.99;
  const static double theta0 = 0.5;
  const static double omega = 1.5;
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
      return sqrt( ((1 - 3*p[0]*p[0]) * p[1] * (1 - p[1]*p[1])) * ((1 - 3*p[0]*p[0]) * p[1] * (1 - p[1]*p[1])) + ((1 - 3*p[1]*p[1]) * p[0] * (1 - p[0]*p[0])) * ((1 - 3*p[1]*p[1]) * p[0] * (1 - p[0]*p[0])) );
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
        values[0] = (1 - 3*p[0]*p[0]) * p[1] * (1 - p[1]*p[1]);
        values[1] = (1 - 3*p[1]*p[1]) * p[0] * (1 - p[0]*p[0]);
      }
      break;
    case 2:
      {
        values[0] = 0.0;
        values[1] = 0.0;
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
        double u_x = (1 - 3*p[0]*p[0]) * p[1] * (1 - p[1]*p[1]);
        double u_y = (1 - 3*p[1]*p[1]) * p[0] * (1 - p[0]*p[0]);
        double u_xx = 6*p[0]*p[1]*( p[1]*p[1] - 1);
        double u_yy = 6*p[1]*p[0]*( p[0]*p[0] - 1);
        double u_xy = (1 - 3*p[0]*p[0]) * (1 - 3*p[1]*p[1]);
        double nabla_u_square = u_x*u_x + u_y*u_y;
        double r;

        r = (param_p-2) * pow(nabla_u_square, (param_p-4)*0.5) * (u_x*u_x*u_xx + u_y*u_y*u_yy + 2*u_x*u_y*u_xy)
          + pow(nabla_u_square, (param_p-2)*0.5) * (u_xx + u_yy);
        return -r;
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
        double u_x = (1 - 3*p[0]*p[0]) * p[1] * (1 - p[1]*p[1]);
        double u_y = (1 - 3*p[1]*p[1]) * p[0] * (1 - p[0]*p[0]);
        double nabla_u_square = u_x*u_x + u_y*u_y;
        double factor = pow(nabla_u_square, (param_p-2)*0.5);

        values[0] = factor * u_x;
        values[1] = factor * u_y;
      }
      break;
      case 2:
      {
        double u_x = d_exact_solution->value(p,0);
        double u_y = d_exact_solution->value(p,1);
        double nabla_u_square = u_x*u_x + u_y*u_y;
        double factor;

        if (nabla_u_square < 1.0e-100)
        {
          values[0] = 0;
          values[1] = 0;
        }
        else
        {
          factor = pow(nabla_u_square, (param_p-2)*0.5);
          values[0] = factor * u_x;
          values[1] = factor * u_y;
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
  Function<2>* d_exact_solution;
};


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



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

  //! Constructor
  Coeff_pPoisson(const AggregatedFrame<IBASIS, DIM, DIM>* frame)
      : frame_(frame)
  {

  }

  double param_p;

  //! point evaluation of ( min{ max{cutoff_small, |\nabla a(x)|}, cutoff_big} )^{p-2}
  double value(const Point<DIM>& p, const unsigned int component = 0) const
  {

    double nabla_square_norm(0.);
    double r;

    // measure time for function evaluation
    clock_t begin_EVALUATE_a = clock();

    //_frame_eval.evaluate_deriv( *frame_, a_coeffs, p, deriv_values);
    _frame_eval.evaluate_deriv_tuned( *frame_, a_coeffs, p, deriv_values);

    clock_t end_EVALUATE_a = clock();
    double elapsed_secs = double(end_EVALUATE_a - begin_EVALUATE_a) / CLOCKS_PER_SEC;
    time_consumption_of_coeff_func += elapsed_secs;


#if 0
//! sanity-check for evaluate

    FixedArray1D<double,DIM> deriv_values2;
    double diff0, diff1;

   _frame_eval.evaluate_deriv( *frame_, a_coeffs, p, deriv_values2);

    diff0 = deriv_values[0] - deriv_values2[0];
    diff1 = deriv_values[1] - deriv_values2[1];

    cout << setprecision(20);
    if (abs(diff0) > 1e-15)
    {
      cout << "ERROR: deriv_values mismatch in x!!" << endl;
      cout << deriv_values[0] << " and " << deriv_values2[0] << endl;
      exit(1);
    }

    if (abs(diff1) > 1e-15)
    {
      cout << "ERROR: deriv_values mismatch in y!!" << endl;
      cout << deriv_values[1] << " and " << deriv_values2[1] << endl;
      exit(1);
    }
    cout << setprecision(6);
//! sanity-check ENDE
#endif

    //! compute l_2-norm of gradient
    for (unsigned int i = 0; i < DIM; i++)
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

  //! set relaxation parameters 'cutoff_small' and 'cutoff_big', and set coefficients 'a_coeffs' of the wavelet frame expansion of a(x)
  void setcoeff(const double epsilon_lb, const double epsilon_ub, const InfiniteVector<double, typename AggregatedFrame<IBASIS, DIM, DIM>::Index> coeffs)
  {
    cutoff_small = epsilon_lb;
    cutoff_big = epsilon_ub;
    a_coeffs.clear();
    a_coeffs += coeffs;
  }

  //! data
  double cutoff_small, cutoff_big;
  mutable FixedArray1D<double,DIM> deriv_values;
  //FixedArray1D<double,DIM> deriv_values2;
  InfiniteVector<double, typename AggregatedFrame<IBASIS, DIM, DIM>::Index> a_coeffs;
  const AggregatedFrame<IBASIS, DIM, DIM>* frame_;
  EvaluateFrame<IBASIS, DIM, DIM> _frame_eval;
};

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



int main()
{
  cout << "Kacanov-type iteration method for the p-Poisson problem on the L-shaped domain ('Mult_Schw' is used for solution of the linear subproblems) ..." << endl;


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                   Parameters and flags
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  const int PARAM_EXAMPLE = 2;               // Encodes the p-Poisson equation to be solved (all on the L-domain with homog. Dirichlet bc's). The corresponding right-hand
                                             // side is defined in class 'TestRHS' (and TestRHS_AR resp.), the explicit solution (if known) in class 'Exact_Solution'.
                                             // '1' corresponds to example 'Bump' (which admits a smooth polynomial solution),
                                             // '2' corresponds to example 'Singular'

  const int DIM = 2;                         // space dimension of the p-Poisson equation

  const double param_p = 1.5;                // the parameter 'p' of the p-Poisson problem

  const int d  = 3;
  const int dT = 3;

  //const int MAX_LEV_WAV_ALL = 7;           // the maximal level of the wavelets
const int MAX_LEV_WAV_ALL = 4; // few levels for testing (speedup)

  const int max_lev_wav_basis = MAX_LEV_WAV_ALL;           // the maximal level of the wavelet bases (CUBEBASIS) which constitute the frame
  const int max_lev_wav_rhs = MAX_LEV_WAV_ALL;             // the maximal level of the wavelets considered by 'compute_rhs()' for precomputation of the rhs coefficients
  const int max_lev_wav_CDD1 = MAX_LEV_WAV_ALL;            // the maximal level of the wavelets considered by 'CDD1_SOLVE'
  const int max_lev_wav_u_exact = MAX_LEV_WAV_ALL;         // the maximal level of the wavelets considered for the expansion of the exact solution u_exact

  //const int max_kacanov_iterations = 15;     // maximal number of Kacanov iterations. in each iteration we solve a lin. elliptic PDE with nonconst. coeff.
const int max_kacanov_iterations = 2; // few iterations for testing (speedup)

  //const int min_res_quad_a = 8;              // parameter for the evaluation of the bilinear form a(,). maximal side length (=2^-min_resolution_quadrature) of the cubes
                                               // used for the composite Gauss quadrature
const int min_res_quad_a = 2; // bad resolution for testing (speedup)

  const int doe_quad_a = 3;                  // degree of exactness of the Gauss quadrature rule used in 'a' (Note: only odd doe possible. For even value, doe will be value-1 !!)


  //const int min_res_quad_f = 9;              // parameter for the evaluation of the functional f(). maximal side length (=2^-min_resolution_quadrature_f) of the cubes
                                               // used for the composite Gauss quadrature. Set value > 0 if rhs is not very smooth.
const int min_res_quad_f = 2; // bad resolution for testing (speedup)

  const int doe_quad_f = 5;                  // degree of exactness of the Gauss quadrature rule used in 'f' (Note: only odd doe possible. For even value, doe will be value-1 !!)


  const int grid_resolution_sm = 6;          // resolution of sampled mappings (and matlab plots)


  const int STRATEGY_EPSILON_STABILIZATION = 1;            // flag for the update strategy of 'epsilon_stabilization': 0 = constant during the whole Kacanov iteration,
                                                           // 1 = algebraic decay, 2 = exponential decay

  const int STRATEGY_EPSILON_MS = 1;                       // flag for the update strategy of 'epsilon_MS': 0 = constant during the whole Kacanov iteration,
                                                           // 1 = exponential decay

  double epsilon_stabilization_lower = 1.0;                                  // initial value of the stabilization parameter for the Kacanov iteration, i.e.
  double epsilon_stabilization_upper = (1.0/epsilon_stabilization_lower);     // $-\div( ( epsilon_stabilization_lower < \lvert \nabla u_n \rvert < epsilon_stabilization_upper )^(p-2) \nabla u_{n+1} ) = f$

 // double epsilon_MS = 0.01;                      // The target $\ell_2$-accuracy for the multiplicative Schwarz algorithm 'Mult_Schw'. Initial value.
double epsilon_MS = 1.01;  // bad accuracy for testing (speedup)

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                          Initialization of external variables
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//! initialize time measurement
  time_consumption_of_a = 0.0;
  time_consumption_of_a_same_patches = 0.0;
  time_consumption_of_a_different_patches = 0.0;
  time_consumption_of_NGROW = 0.0;
  time_consumption_of_GALERKIN = 0.0;
  time_consumption_of_NRESIDUAL = 0.0;
  time_consumption_of_coeff_func = 0.0;
  time_consumption_of_calc_deriv = 0.0;
  time_consumption_of_calc_deriv_outer = 0.0;
  time_consumption_of_compute_rhs = 0.0;
  time_consumption_of_compute_diagonal = 0.0;
  time_consumption_of_eval_u_eps = 0.0;
  time_consumption_of_eval_deriv_u_eps = 0.0;

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                     typedef's
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  typedef PBasis<d,dT> Basis1D;             // Primbs basis on the intervall [0,1]
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
//  typedef CubeBasis<Basis1D,2> Basis;       // tensor product wavelet basis on [0,1]^2
//  typedef MappedCubeBasis<Basis1D,2,2> MappedBasis;
  typedef Frame2D::Index Index;
//  typedef CubeIndex<Basis1D,2,MappedCubeBasis<Basis1D,2,2> > CIndex;

//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                 Variable declarations
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//! -------------------------------------------------------------------------------------------------------------
//! preparations for the definition of the frame (charts, atlas, bc's)


  //! -------------------------------------     the charts     --------------------------------------------------

  //Matrix<double> A(DIM,DIM); // for AffineLinearMapping
  //A(0,0) = 1. + OVERLAP;
  //A(1,1) = 1.0;
  FixedArray1D<double,DIM> A;  // for SimpleAffineLinearMapping
  A[0] = 1. + OVERLAP;
  A[1] = 1.0;
  Point<2> b;
  b[0] = -OVERLAP;
  b[1] = -1.;
  SimpleAffineLinearMapping<2> affineP(A,b);

  //Matrix<double> A2(DIM,DIM);  // for AffineLinearMapping
  //A2(0,0) = 1.0;
  //A2(1,1) = 2.0;
  FixedArray1D<double,DIM> A2;  // for SimpleAffineLinearMapping
  A2[0] = 1.0;
  A2[1] = 2.0;
  Point<2> b2;
  b2[0] = -1.0;
  b2[1] = -1.0;
  SimpleAffineLinearMapping<2> affineP2(A2,b2);

  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
  charts[1] = &affineP2;

  // adjacency matrix of the patches
  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;

  //! --------------------------------------     the atlas     ---------------------------------------------------

  Atlas<DIM,DIM> Lshaped(charts,adj);
  cout << Lshaped << endl;

  //! --------------------------------------     the bc's     ----------------------------------------------------

  Array1D<FixedArray1D<int,2*DIM> > bc(2);      // to specify primal boundary conditions

  //! primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 1;

  bc[0] = bound_1;

  //! primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;
  bound_2[2] = 1;
  bound_2[3] = 1;

  bc[1] = bound_2;

  Array1D<FixedArray1D<int,2*DIM> > bcT(2);     //to specify dual boundary conditions

  //! dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;
  bound_3[2] = 0;
  bound_3[3] = 0;

  bcT[0] = bound_3;

  //! dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
  bound_4[2] = 0;
  bound_4[3] = 0;

  bcT[1] = bound_4;


  //! -----------------------------------------     the frame     -----------------------------------------------------

  //AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT, 6);
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, max_lev_wav_basis);


  //! ------------------------------------     the rhs & exact solution     -------------------------------------------

  // constant rhs
  Vector<double> value(1);
  value[0] = 1;
  ConstantFunction<DIM> const_fun(value);


#if EXAMPLE_SINGULARITY
  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;

  //! exact solution
  CornerSingularity u_exact(origin, 0.5, 1.5);
  CornerSingularityGradient deriv_u_exact(origin, 0.5, 1.5);

  //! right-hand side f_ar: \Omega -> \R^DIM (alternative representation of the rhs functional),
  //! i.e., the wavelet coefficients of f take the form: rhs_{\lambda} = \int_{\Omega} f_ar(x) * \nabla Psi_{\lambda}(x) dx
  TestRHS_AR<PARAM_EXAMPLE> rhs_ar;
  rhs_ar.param_p = param_p;
  rhs_ar.d_exact_solution = &deriv_u_exact;

//! right-hand side f:
  TestRHS_AR<PARAM_EXAMPLE> rhs;    // DUMMY !! Use only with #define PPE_USE_ALTERNATIVE_F 1  AND  #define PPDF_USE_ALTERNATIVE_F 1  !!!
  rhs.param_p = param_p;

#else
  //! right-hand side f: \Omega -> \R (as pointwise function),
  //! i.e.: rhs_{\lambda} = \int_{\Omega} f(x) * Psi_{\lambda}(x) dx
  TestRHS<PARAM_EXAMPLE> rhs;
  rhs.param_p = param_p;

  //! right-hand side f_ar: \Omega -> \R^d (alternative representation of the rhs functional),
  //! i.e., the wavelet coefficients of f take the form: rhs_{\lambda} = \int_{\Omega} f_ar(x) * \nabla Psi_{\lambda}(x) dx
  TestRHS_AR<PARAM_EXAMPLE> rhs_ar;
  rhs_ar.param_p = param_p;

  //! exact solution
  Exact_Solution<PARAM_EXAMPLE> u_exact;                            // exact solution as a function
  Deriv_Exact_Solution<PARAM_EXAMPLE> deriv_u_exact;                // derivative of the exact solution as a function
#endif

  //! ------------------------------     the coefficient of the p-Laplace operator     --------------------------------------

  Coeff_pPoisson<Basis1D,2> mein_coeff_pP(&frame);
  mein_coeff_pP.param_p = param_p;
  mein_coeff_pP.setcoeff(epsilon_stabilization_lower, epsilon_stabilization_upper, InfiniteVector<double, Index>());


  //! ---------------------------     the continuous bvp (with nonconst. coefficients)     ----------------------------------

  PoissonBVP_Coeff<2> poisson_coeff(&mein_coeff_pP, &rhs, &rhs_ar);
  // PoissonBVP_Coeff<2> poisson_coeff(&mein_coeff_pP, &rhs);
  // PoissonBVP_Coeff<2> poisson_coeff(&mein_coeff_pP, &const_fun);  // rhs = const


  //! --------------------------------------     the discrete problem     -------------------------------------------

  pPoissonEquation<Basis1D,DIM> discrete_poisson(&poisson_coeff, &frame, param_p,
                                               max_lev_wav_rhs,
                                               doe_quad_a, min_res_quad_a,
                                               doe_quad_f, min_res_quad_f);


  //! ---------------------------------------     the cached problem     --------------------------------------------

#if CP_USE_PRECOMPUTED_NORM

  //  L-shaped: (-1,1)x(-1,0) \cup (-1,0)x(-1,1), PBasis

  //   // (d,dT) = (2,2)
  //   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.0225, 1.0/0.146);
  //   discrete_poisson.set_norm_A(5.0225);
  //   // optimistic guess:
  //   discrete_poisson.set_Ainv(1.0/0.146);

  // (d,dT) = (3,5)
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  //CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  //CachedProblemLocal<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  //CachedProblemLocal<pPoissonEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);

  // (d,dt) = (3,3), p=1.5, example 1 ('bump')
  CachedProblemLocal<pPoissonEquation<Basis1D,DIM> > problem(&discrete_poisson, 7.3, 1.0/0.146);

  problem.c1_patch0 = 0.05;
  problem.c2_patch0 = 4.9;
  problem.c1_patch1 = 0.05;
  problem.c2_patch1 = 4.9;


#else
  CachedProblemLocal<pPoissonEquation<Basis1D,DIM> > problem(&discrete_poisson);
  problem.set_normAinv(1.0/0.146); // quick hack
#endif


  //! -------------------------------------     computation of constants     --------------------------------------------
  // estimate extremal eigenvalues of local stiffness matrices and largest eigenvalue of whole stiffness matrix

#if COMPUTECONSTANTS

  set<Index> Lambda_0;
  set<Index> Lambda_1;
  set<Index> Lambda;
  for (FrameIndex<Basis1D,DIM,DIM> lambda = FrameTL::first_generator<Basis1D,DIM,DIM,Frame2D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,DIM,DIM,Frame2D>(&frame, max_lev_wav_basis); ++lambda)
  {
    Lambda.insert(lambda);
    if (lambda.p() == 0)
    {
      Lambda_0.insert(lambda);
    }
    else
    {
      Lambda_1.insert(lambda);
    }
  }

  SparseMatrix<double> stiff;

  // starting vector for Power and Inverse Power Iteration
  Vector<double> x(Lambda_0.size());
  x = 1;
  // number of iterations in Power and Inverse Power Iteration
  unsigned int iter= 0;

  WaveletTL::setup_stiffness_matrix(problem, Lambda_0, stiff);

  cout << "computing smallest eigenvalue of stiffness matrix on patch 0" << endl;
  double lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "smallest eigenvalue of stiffness matrix on patch 0 is " << lmin << endl;
  cout << "iterations:  " << iter << endl;

  cout << "computing largest eigenvalue of stiffness matrix on patch 0" << endl;
  double lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of stiffness matrix on patch 0 is " << lmax << endl;
  cout << "iterations:  " << iter << endl;

  WaveletTL::setup_stiffness_matrix(problem, Lambda_1, stiff);

  x.resize(Lambda_1.size()); x = 1;
  cout << "computing smallest eigenvalue of stiffness matrix on patch 1" << endl;
  lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "smallest eigenvalue of stiffness matrix on patch 1 is " << lmin << endl;
  cout << "iterations:  " << iter << endl;

  cout << "computing largest eigenvalue of stiffness matrix on patch 1" << endl;
  lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of stiffness matrix on patch 1 is " << lmax << endl;
  cout << "iterations:  " << iter << endl;

  x.resize(Lambda.size()); x = 1;
  WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);

  cout << "computing smallest eigenvalue of whole stiffness matrix" << endl;
  lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "smallest eigenvalue of whole stiffness matrix is " << lmin << endl;
  cout << "iterations:  " << iter << endl;

  cout << "computing largest eigenvalue of whole stiffness matrix" << endl;
  lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of whole stiffness matrix is " << lmax << endl;
  cout << "iterations:  " << iter << endl;


  // (d,dt) = (2,2), jmax = 5, jmin = 3:
  // patch 0: \lambda_{\min}^0 = 0.0961009, \lambda_{\max}^0 = 3.17632
  // patch 1: \lambda_{\min}^1 = 0.0961009, \lambda_{\max}^1 = 3.17632
  // whole domain: \lambda_{\max} = 5.01773

  // (d,dt) = (2,2), jmax = 5, jmin = 4:
  // patch 0: \lambda_{\min}^0 = 0.0277685, \lambda_{\max}^0 = 3.02486
  // patch 1: \lambda_{\min}^1 = 0.0277685, \lambda_{\max}^1 = 3.02486
  // whole domain: \lambda_{\max} = 4.60975

  // (d,dt) = (3,3), jmax = 5, jmin = 3:
  // patch 0: \lambda_{\min}^0 = 0.0748624, \lambda_{\max}^0 = 4.74753
  // patch 1: \lambda_{\min}^1 = 0.0748624, \lambda_{\max}^1 = 4.74753
  // whole domain: \lambda_{\max} = 6.98681

  // (d,dt) = (3,3), jmax = 5, jmin = 4:
  // patch 0: \lambda_{\min}^0 = 0.0664664, \lambda_{\max}^0 = 4.76616
  // patch 1: \lambda_{\min}^1 = 0.0664664, \lambda_{\max}^1 = 4.76616
  // whole domain: \lambda_{\max} = 6.98986

    // (d,dt) = (4,4), jmax = 5, jmin = 4:
  // patch 0: \lambda_{\min}^0 = 0.00285239 , \lambda_{\max}^0 = 8.5811
  // patch 1: \lambda_{\min}^1 = 0.00285239  \lambda_{\max}^1 = 8.57653
  // whole domain: \lambda_{\max} = 12.3335

  abort();

#endif


  //! -------------------------------------     data class     --------------------------------------------
  // for computation and storage of: the approximation error (in various norms), plots, parameters, energy etc.

  pPoissonDataFrame< Basis1D, CachedProblemLocal<pPoissonEquation<Basis1D,DIM> > > p_poisson_data(&frame, &problem, &u_exact, &deriv_u_exact, grid_resolution_sm);


  //! compute the energy of the exact solution and save a (matlab) function plot
  clock_t begin_precomp_exact_solution = clock();

  p_poisson_data.precompute_exact_solution(5, 8);    // example bump: since exact solution is smooth (polynomial), we use (doe, min_res_quad) = (5,8)

  clock_t end_precomp_exact_solution = clock();
  double time_consumption_of_precomp_exact_solution = double(end_precomp_exact_solution - begin_precomp_exact_solution) / CLOCKS_PER_SEC;
  cout << "time for precompute_exact_solution: " << time_consumption_of_precomp_exact_solution << endl;

  //! -----------------------------------     the approximants     ------------------------------------------

  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);       // stores current local approximations on each patch, as well as the current global approximation

  InfiniteVector<double, Index> u_epsilon_old;                  // global approximation of last iteration (as a guess for the next iterate)
  u_epsilon_old.clear();


//! +++++++++++++++++++++++++++++++++++++++++++++++++++++     *END* of Variable declarations     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                  Start of the Kacanov iteration for solving p-Poisson
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#if 0
//! parameter test (doe, min_res) for a:
double dummy_a1 = problem.a(*(discrete_poisson.frame().get_wavelet(1)),*(discrete_poisson.frame().get_wavelet(8)));
double dummy_a2 = problem.a(*(discrete_poisson.frame().get_wavelet(99)),*(discrete_poisson.frame().get_wavelet(100)));
double dummy_a3 = problem.a(*(discrete_poisson.frame().get_wavelet(400)),*(discrete_poisson.frame().get_wavelet(405)));
double dummy_a4 = problem.a(*(discrete_poisson.frame().get_wavelet(1000)),*(discrete_poisson.frame().get_wavelet(1004)));
double dummy_a5 = problem.a(*(discrete_poisson.frame().get_wavelet(1900)),*(discrete_poisson.frame().get_wavelet(1915)));

cout << "(1,8) : " << *(discrete_poisson.frame().get_wavelet(1)) << " , " << *(discrete_poisson.frame().get_wavelet(8)) << endl;
cout << "(99,100) : " << *(discrete_poisson.frame().get_wavelet(99)) << " , " << *(discrete_poisson.frame().get_wavelet(100)) << endl;
cout << "(400,405) : " << *(discrete_poisson.frame().get_wavelet(400)) << " , " << *(discrete_poisson.frame().get_wavelet(405)) << endl;
cout << "(1000,1004) : " << *(discrete_poisson.frame().get_wavelet(1000)) << " , " << *(discrete_poisson.frame().get_wavelet(1004)) << endl;
cout << "(1900,1915) : " << *(discrete_poisson.frame().get_wavelet(1900)) << " , " << *(discrete_poisson.frame().get_wavelet(1915)) << endl << endl;

cout << std::setprecision(12);

cout << "a(1,8) = " << dummy_a1 << endl;
cout << "a(99,100) = " << dummy_a2 << endl;
cout << "a(400,405) = " << dummy_a3 << endl;
cout << "a(1000,1004) = " << dummy_a4 << endl;
cout << "a(1900,1915) = " << dummy_a5 << endl;

cout << std::setprecision(6);

exit(1);
#endif // end parameter test

  for (int i = 0; i < max_kacanov_iterations; i++)
  {


  //! ----------------------------------     Approximate min/max eigenvalues of stiffness matrix (local and globaL)     -------------------------------------------
  // determine c1,c2 (i.e. min., max. eigenvalue) of local stiffness matrix for each patch, as well as norm of whole stiffness matrix A

#if !CP_USE_PRECOMPUTED_NORM
    set<Index> Lambda_0;
    set<Index> Lambda_1;
    set<Index> Lambda;
    for (Index lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&(problem.basis()), problem.basis().j0());
         lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&(problem.basis()), problem.basis().jmax()); ++lambda)
    {
      Lambda.insert(lambda);
      if (lambda.p() == 0)
      {
        Lambda_0.insert(lambda);
      }
      else
      {
        Lambda_1.insert(lambda);
      }
    }

    SparseMatrix<double> stiff;

    // starting vector for Power and Inverse Power Iteration
    Vector<double> x(Lambda_0.size());
    x = 1;
    // number of iterations in Power and Inverse Power Iteration
    unsigned int iter= 0;

    WaveletTL::setup_stiffness_matrix(problem, Lambda_0, stiff);

    cout << "computing smallest eigenvalue of stiffness matrix on patch 0" << endl;
    double lmin_0 = InversePowerIteration(stiff, x, 1e-6, 1000, iter);
    cout << "smallest eigenvalue of stiffness matrix on patch 0 is " << lmin_0 << endl;
    cout << "iterations:  " << iter << endl;

    cout << "computing largest eigenvalue of stiffness matrix on patch 0" << endl;
    double lmax_0 = PowerIteration(stiff, x, 1e-6, 1000, iter);
    cout << "largest eigenvalue of stiffness matrix on patch 0 is " << lmax_0 << endl;
    cout << "iterations:  " << iter << endl;

    WaveletTL::setup_stiffness_matrix(problem, Lambda_1, stiff);

    x.resize(Lambda_1.size()); x = 1;
    cout << "computing smallest eigenvalue of stiffness matrix on patch 1" << endl;
    double lmin_1 = InversePowerIteration(stiff, x, 1e-6, 1000, iter);
    cout << "smallest eigenvalue of stiffness matrix on patch 1 is " << lmin_1 << endl;
    cout << "iterations:  " << iter << endl;

    cout << "computing largest eigenvalue of stiffness matrix on patch 1" << endl;
    double lmax_1 = PowerIteration(stiff, x, 1e-6, 1000, iter);
    cout << "largest eigenvalue of stiffness matrix on patch 1 is " << lmax_1 << endl;
    cout << "iterations:  " << iter << endl;

    x.resize(Lambda.size()); x = 1;
    WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);

    cout << "computing smallest eigenvalue of whole stiffness matrix" << endl;
    double lmin = InversePowerIteration(stiff, x, 1e-6, 1000, iter);
    cout << "smallest eigenvalue of whole stiffness matrix is " << lmin << endl;
    cout << "iterations:  " << iter << endl;

    cout << "computing largest eigenvalue of whole stiffness matrix" << endl;
    double lmax = PowerIteration(stiff, x, 1e-6, 1000, iter);
    cout << "largest eigenvalue of whole stiffness matrix is " << lmax << endl;
    cout << "iterations:  " << iter << endl;


    cout << endl << "MultSchw parameters stiffnes matrix:" << endl
         << "  lmin_0 = " << lmin_0 << endl
	     << "  lmax_0 = " << lmax_0 << endl
	     << "  lmin_1 = " << lmin_1 << endl
	     << "  lmax_1 = " << lmax_1 << endl
	     << "  lmin = " << lmin << endl
	     << "  lmax = " << lmax << endl << endl;


    problem.c1_patch0 = lmin_0;
    problem.c2_patch0 = lmax_0;
    problem.c1_patch1 = lmin_1;
    problem.c2_patch1 = lmax_1;

    problem.set_normA(lmax);

#endif




  //! --------------------------------------     Compute next iterate u_epsilon     ------------------------------------------------


    cout << "adaptive multiplicative Schwarz started..." << endl;

    cout << "p_Poisson: i=" << i << ", Stabilization parameter: epsilon_lower =  " << epsilon_stabilization_lower << ", epsilon_upper =  " << epsilon_stabilization_upper << endl;

    clock_t begin_Mult_Schw = clock();

    //! the adaptive multiplicative Schwarz solver
    MultSchw(problem, epsilon_MS, approximations, u_epsilon_old);
    //MultSchw(problem, epsilon_MS, approximations, InfiniteVector<double, Index>());    // to estimate rho: start with zero-guess

    clock_t end_Mult_Schw = clock();
    double time_consumption_of_Mult_Schw = double(end_Mult_Schw - begin_Mult_Schw) / CLOCKS_PER_SEC;
    cout << "adaptive multiplicative Schwarz done." << endl;


    //! store current global approximation for the initial guess in the next iteration
    u_epsilon_old.clear();
    u_epsilon_old += approximations[frame.n_p()];


    //! scale approximations
    for (int j = 0; j <= frame.n_p(); j++)
    {
      approximations[j].scale(&discrete_poisson,-1);
    }


  //! ------------------------------------------     Output time measurement     -----------------------------------------------------

    cout << endl;
    cout << "time for compute_rhs, i = " << i << " : " << time_consumption_of_compute_rhs << endl;
    cout << "time for compute_diagonal, i = " << i << " : " << time_consumption_of_compute_diagonal << endl;
    cout << "time for Mult_Schw, i = " << i << " : " << time_consumption_of_Mult_Schw << endl;
    cout << "time for NGROW, i = " << i << " : " << time_consumption_of_NGROW << endl;
    cout << "time for NRESIDUAL, i = " << i << " : " << time_consumption_of_NRESIDUAL << endl;
    cout << "time for GALERKIN, i = " << i << " : " << time_consumption_of_GALERKIN << endl;
    cout << "time for a() (cached), i = " << i << " : " << time_consumption_of_a << endl;
    cout << "time for a(), same patches, i = " << i << " : " << time_consumption_of_a_same_patches << endl;
    cout << "time for a(), different patches, i = " << i << " : " << time_consumption_of_a_different_patches << endl;
    cout << "time for coeff_func, i = " << i << " : " << time_consumption_of_coeff_func << endl;
    cout << "time for calc_deriv + outer loop, i = " << i << " : " << time_consumption_of_calc_deriv_outer << endl;
    cout << "time for calc_deriv, i = " << i << " : " << time_consumption_of_calc_deriv << endl;
    cout << "time for eval_u_eps, i = " << i << " : " << time_consumption_of_eval_u_eps << endl;
    cout << "time for eval_deriv_u_eps, i = " << i << " : " << time_consumption_of_eval_deriv_u_eps << endl;


  //! -------------------     Save computed data (approximations, parameters etc.) and compute and save errors     -------------------------


    clock_t begin_compute_data = clock();

    cout << "'p_poisson_data': compute_data .." << endl;

    p_poisson_data.compute_data(approximations, param_p, 5, 8, false, epsilon_stabilization_lower, epsilon_stabilization_upper, 0);   // example bump: since exact solution is smooth (polynomial), we use (doe, min_res_quad) = (5,8)

    cout << "'p_poisson_data': .. finished!" << endl;

    clock_t end_compute_data = clock();
    double time_consumption_of_compute_data = double(end_compute_data - begin_compute_data) / CLOCKS_PER_SEC;

    cout << "time for eval_u_eps, i = " << i << " : " << time_consumption_of_eval_u_eps << endl;
    cout << "time for eval_deriv_u_eps, i = " << i << " : " << time_consumption_of_eval_deriv_u_eps << endl;

    p_poisson_data.Data[p_poisson_data.iterations() - 1].epsilon_MS = epsilon_MS;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].epsilon_stabilization = epsilon_stabilization_lower;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].max_lev_wav_CDD1 = max_lev_wav_CDD1;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].max_lev_wav_rhs = max_lev_wav_rhs;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].max_lev_wav_u_exact = max_lev_wav_u_exact;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].min_res_quad_a = min_res_quad_a;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].doe_quad_a = doe_quad_a;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].min_res_quad_f = min_res_quad_f;
    p_poisson_data.Data[p_poisson_data.iterations() - 1].doe_quad_f = doe_quad_f;
   // p_poisson_data.Data[p_poisson_data.iterations() - 1].gamma_CDD1 = gamma_CDD1;

    cout << "time for compute_data, i = " << i << " : " << time_consumption_of_compute_data << endl;



  //! --------------------------------------      Perform several updates before next iteration     -----------------------------------------------


    //! initialize time measurement
    time_consumption_of_a = 0.0;
    time_consumption_of_a_same_patches = 0.0;
    time_consumption_of_a_different_patches = 0.0;
    time_consumption_of_NGROW = 0.0;
    time_consumption_of_GALERKIN = 0.0;
    time_consumption_of_NRESIDUAL = 0.0;
    time_consumption_of_coeff_func = 0.0;
    time_consumption_of_calc_deriv = 0.0;
    time_consumption_of_calc_deriv_outer = 0.0;
    time_consumption_of_compute_rhs = 0.0;
    time_consumption_of_compute_diagonal = 0.0;
    time_consumption_of_eval_u_eps = 0.0;
    time_consumption_of_eval_deriv_u_eps = 0.0;

    if (i == (max_kacanov_iterations-1))
    {
      break;
    }

    //! set minimal resolution of composite gauss quadrature to: max{ (highest wavelet level of u_epsilon) + 1, min_res_quadrature_a }
    std::set<Index> supp_u;
    approximations[frame.n_p()].support(supp_u);
    std::set<Index>::const_iterator it = supp_u.end();

    if (supp_u.empty())
    {
      if (min_res_quad_a > frame.j0())
      {
        discrete_poisson.min_res_quadrature_a = min_res_quad_a;
      }
      else
      {
        discrete_poisson.min_res_quadrature_a = frame.j0() + 1;
      }
    }
    else
    {
      --it;
      Index dummy_index = *it;

      if (min_res_quad_a > dummy_index.j())
      {
        discrete_poisson.min_res_quadrature_a = min_res_quad_a;
      }
      else
      {
        discrete_poisson.min_res_quadrature_a = dummy_index.j() + 1;
      }

      cout << "p_Poisson: last = " << dummy_index << " level = " << dummy_index.j() << endl;
    }


    //! some output..
    cout << "p_Poisson: i=" << i << ", min_res_quadrature_a: " << discrete_poisson.min_res_quadrature_a << endl;


    //! update stabilization parameter
    if (STRATEGY_EPSILON_STABILIZATION == 0)
    {
      ;
    }
    else if (STRATEGY_EPSILON_STABILIZATION == 1)
    {
      epsilon_stabilization_lower = pow(i+2,-1.0);
      epsilon_stabilization_upper = pow(i+2, 1.0);
    }
    else if (STRATEGY_EPSILON_STABILIZATION == 2)
    {
      epsilon_stabilization_lower *= 0.5;
      epsilon_stabilization_upper *= 2.0;
    }


    //! update MultSchwarz accuracy
    if (STRATEGY_EPSILON_MS == 0)
    {
      ;
    }
    else if (STRATEGY_EPSILON_MS == 1)
    {
      epsilon_MS *= 0.5;
    }


    //! update coefficient a
    mein_coeff_pP.setcoeff(epsilon_stabilization_lower, epsilon_stabilization_upper, approximations[frame.n_p()]);    // update coefficient function
    poisson_coeff.set_a(&mein_coeff_pP);             // upadate bvp


    #if !CP_USE_PRECOMPUTED_NORM
    problem.set_normA(0.);
    problem.set_normAinv(0.);
    problem.set_normAinv(1.0/0.146);
    #endif


    problem.clear_cache();                    // clear the entries cache for A

    discrete_poisson.clear_coeff_cache();               // clear cache for the coefficient function

    discrete_poisson.set_bvp(&poisson_coeff);           // set bvp (-> invokes call of compute_rhs and compute_diagonal)

  }


//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++     *END* of Kacanov iteration     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//!                                                                    Save Data / Output
//! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  cout << "Basis1D: primal_polynomial_degree: " << Basis1D::primal_polynomial_degree() << endl;


  //! ----------------------------      save all parameters, computed errors etc. in file. Output to terminal as well     ----------------------------------


  std::ofstream ofs_pPoisson_data("./results/pPoisson_data.txt");                //! open file

  ofs_pPoisson_data << "p-Poisson equation: Results of example " << PARAM_EXAMPLE << " with p = " << param_p << endl << endl;
  for (int i = 0; i < p_poisson_data.iterations(); i++)
  {
    p_poisson_data.print_data(i, std::cout);                                     //! print to terminal
    p_poisson_data.print_data(i, ofs_pPoisson_data);                             //! print to file
  }
  ofs_pPoisson_data.close();                                                     //! close file


  //! ---------------------------------------------      save some error plots (matlab)     ------------------------------------------------


  std::ofstream ofs_plot_error_energy("./results/error_energy.m");             //! open file
  std::ofstream ofs_plot_error_W1p("./results/error_W1p.m");                   //! open file
  std::ofstream ofs_plot_error_W1_infty("./results/error_W1_infty.m");         //! open file

  p_poisson_data.save_error_plots(ofs_plot_error_energy,                       //! save error plots
                                  ofs_plot_error_W1p,
                                  ofs_plot_error_W1_infty);

  ofs_plot_error_energy.close();                                               //! close file
  ofs_plot_error_W1p.close();                                                  //! close file
  ofs_plot_error_W1_infty.close();                                             //! close file


  //! -----------------------------------------------      save plot of rhs (matlab)     ---------------------------------------------------


  //! rhs as sampled mapping
  Array1D<SampledMapping<2> > rhs_sm(frame.n_p());

  for (int i = 0; i < frame.n_p(); i++)
  {
    rhs_sm[i] = SampledMapping<2>(*(frame.atlas()->charts()[i]), grid_resolution_sm);  // all zero
    rhs_sm[i].add(1.,SampledMapping<2>(rhs_sm[i], rhs ));
  }

  //! rhs as Matlab plot
  std::ofstream rhs_stream("./results/plots/rhs.m");
  matlab_output(rhs_stream, rhs_sm);
  rhs_stream.close();


//! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     *END* of Save Data / Output     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  return 0;

}  //! *END* of main()



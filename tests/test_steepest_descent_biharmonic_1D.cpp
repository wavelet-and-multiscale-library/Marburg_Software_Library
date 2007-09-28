//#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
//#define _WAVELETTL_CACHEDPROBLEM_VERBOSITY 2
//#define _FRAMETL_SIMPLEBIHARMONICEQUATION_VERBOSITY 2
//#define _WAVELETTL_INDEX_VERBOSITY 1

// choose which righthand side to use
//#define CONSTANT_RHS

// choose which basis to use
#define BASIS_P

//#ifdef BASIS_S
//#define MAX_LOOPS 10000
//#else
//#define MAX_LOOPS 7000
//#endif

// tweaking the adaptive algorithm
#ifdef BASIS_S
 #define ETA_STEP 0.9995;
#else
 #define ETA_STEP 0.995;
#endif

#define SAVE_ASYMPTOTIC 10
#define SAVE_LOG

#include <fstream>
#include <iostream>
#include <time.h> 

#ifdef BASIS_S
#include <interval/s_basis.h>
#include <interval/s_support.h>
#include <interval/interval_evaluate.h>
#include <interval/adapted_basis.h>
#include <interval/adapted_support.h>
#define BASIS_NAME "s"
#else
#include <interval/p_basis.h>
#include <interval/p_support.h>
#include <interval/p_evaluate.h>
#define BASIS_NAME "p"
#endif
#include <aggregated_frame.h>
#include <frame_evaluate.h>
#include <frame_support.h>
#include <frame_index.h>

#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <simple_biharmonic_equation.h>
#include <biharmonic_1d_testcase.h>

#include <galerkin/cached_problem.h>

// typedefs for choice of basis
#ifdef BASIS_S
 #define D_PRIMAL 4
 #define D_DUAL 2
 typedef AdaptedBasis<SBasis> Basis1D;
#else
 #define D_PRIMAL 4
 #define D_DUAL 4
 typedef PBasis<D_PRIMAL,D_DUAL> Basis1D;
#endif
typedef CubeBasis<Basis1D,1> IntervalBasis;
typedef AggregatedFrame<Basis1D,1,1> Frame1D;
typedef Frame1D::Index Index;

#include <simplified_steepest_descent.h>

using std::cout;
using std::endl;

using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
 
  cout << "Testing steepest descent with the biharmonic equation in 1D ..." << endl;
  
  const int DIM = 1;
  const int jmax = 13;
  const double epsilon = 1e-3;

  // set up frame ***********************************************************
  // first patch
  Matrix<double> A1(DIM,DIM);
  A1(0,0) = 0.7;
  Point<1> b1;
  b1[0] = 0.;
  AffineLinearMapping<1> affineP1(A1,b1);

  // second patch
  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 0.7;
  Point<1> b2;
  b2[0] = 1-A2.get_entry(0,0);
  AffineLinearMapping<1> affineP2(A2,b2);

  // put things together
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP1;
  charts[1] = &affineP2;

  // adjacency matrix
  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;

  // specify primal boundary the conditions ---------------------------------
  Array1D<FixedArray1D<int,2*DIM> > bc(2);
  
  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 2;
  #ifdef BASIS_S
  bound_1[1] = 2;
  #else
  bound_1[1] = D_PRIMAL - 1;
  #endif
  
  bc[0] = bound_1;
  
  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  #ifdef BASIS_S
  bound_2[0] = 2;
  #else
  bound_2[0] = D_PRIMAL - 1;
  #endif
  bound_2[1] = 2;
  
  bc[1] = bound_2;

  // specify dual boundary the conditions -----------------------------------
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;

  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
 
  bcT[1] = bound_4;

  // atlas with charts ------------------------------------------------------
  Atlas<DIM,DIM> interval(charts,adj);
  cout << interval << endl;

  // construct frame --------------------------------------------------------
  Frame1D frame(&interval, bc, jmax);


  // setup biharmonic equation **********************************************
  #ifdef CONSTANT_RHS
  // constant right-hand sinde, exact solution is 16 x^2 (1-x)^2
  ConstantFunction<DIM> const_fun(Vector<double>(1, "384"));
  Functional<Basis1D,1,1> rhs(&const_fun, &frame); // right-hand side
  Polynomial<double> exact_solution(Vector<double>(5, "0 0 16 -32 16")); // exact solution
  #else
  Biharmonic1D_Solution exact_solution; // exact solution
  Biharmonic1D_RHS<Basis1D> rhs(&frame); // right-hand side
  #endif

  SimpleBiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&rhs, &frame, jmax, TrivialAffine);
 
  CachedProblem<SimpleBiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 5, 1.0/0.146);


  // start adapted solver ***************************************************
  clock_t tstart, tend;
  double time;
  tstart = clock();

  InfiniteVector<double, Index> u_epsilon;
  simplified_steepest_descent_SOLVE(problem, epsilon, u_epsilon, jmax);

  tend = clock();
  time = (double)((tend-tstart)/CLOCKS_PER_SEC);
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "steepest descent done" << endl;

  u_epsilon.scale(&discrete_biharmonic,-1);


  // evaluate solution and do some plotting *********************************
  cout << "Evaluating solution ..." << endl;
  EvaluateFrame<Basis1D,1,1> evalObj;

  Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, u_epsilon, true, 11); // expand in primal basis
  cout << "... finished plotting approximate solution" << endl;
  Array1D<SampledMapping<1> > error = evalObj.evaluate_difference(frame, u_epsilon, exact_solution, 11);
  cout << "... finished plotting error" << endl;
  
  std::ofstream ofs_approx("steep_1D_bi_approx_sol_out.m");
  matlab_output(ofs_approx,U);
  ofs_approx.close();

  std::ofstream ofs_error("steep_1D_bi_error_out.m");
  matlab_output(ofs_error,error);
  ofs_error.close();

  return 0;
}

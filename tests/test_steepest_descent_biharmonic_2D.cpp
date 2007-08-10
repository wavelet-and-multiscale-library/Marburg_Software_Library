#define BASIS_S
#define RHS_QUADRATURE_GRANULARITY 5
#define MAX_LOOPS 12000
#define PLOT_RESOLUTION 6
#define SAVE_ASYMPTOTIC
//#define SAVE_LOG

#ifdef BASIS_S
 #define ETA_STEP 0.995
#else
 #define ETA_STEP 0.99
#endif

#include <fstream>
#include <iostream>
#include <time.h>

#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>

#ifdef BASIS_DS
#include <interval/ds_basis.h>
#define BASIS_NAME "ds"
#endif
#ifdef BASIS_P
#include <interval/p_basis.h>
#include <interval/p_support.h>
#include <interval/p_evaluate.h>
#define BASIS_NAME "p"
#endif
#ifdef BASIS_S
#include <interval/s_basis.h>
#include <interval/s_support.h>
#include <interval/interval_evaluate.h>
#include <interval/adapted_basis.h>
#include <interval/adapted_support.h>
#define BASIS_NAME "s"
#endif

#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <aggregated_frame.h>
#include <frame_support.h>
#include <frame_index.h>
#include <frame_evaluate.h>

#include <numerics/corner_singularity_biharmonic.h>
#include <simple_biharmonic_equation.h>

#include <galerkin/galerkin_utils.h>
#include <galerkin/cached_problem.h>
#include <simplified_steepest_descent.h>

#include <map>

// typedefs for choice of basis
#ifdef BASIS_DS
typedef DSBasis<4,6> Basis1D;
#endif
#ifdef BASIS_P
#define D_PRIMAL 4
#define D_DUAL 4
typedef PBasis<D_PRIMAL,D_DUAL> Basis1D;
#endif
#ifdef BASIS_S
#define D_PRIMAL 4
#define D_DUAL 2
typedef AdaptedBasis<SBasis> Basis1D;
#endif
typedef AggregatedFrame<Basis1D,2,2> Frame2D;
typedef CubeBasis<Basis1D> Basis;
typedef Frame2D::Index Index;


using std::cout;
using std::endl;

using MathTL::ConstantFunction;
using MathTL::CornerSingularityBiharmonicRHS;
using MathTL::CornerSingularityBiharmonic;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;
using WaveletTL::CachedProblem;

using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


/* functions for saving the and loading the (pre)computed right-hand side */
template <class C>
void Rhs_output (const Vector<C> v, const int jmax) 
{
  SparseMatrix<C> S = SparseMatrix<C>(v.size());

  for(unsigned int i=0; i< v.size(); i++)
    S.set_entry( 0, i, v[i]);

  ostringstream filename;
  filename << "Rhs_2D_" << BASIS_NAME << "_jmax" << jmax << "_N" << RHS_QUADRATURE_GRANULARITY << "_out";
  S.matlab_output(filename.str().c_str(), "Matrix", 1);
 }
  
template <class C>
void Rhs_input(Vector<C> &v, const int jmax)
{
  SparseMatrix<C> S;
  ostringstream filename;
  filename << "Rhs_2D_" << BASIS_NAME << "_jmax" << jmax << "_N" << RHS_QUADRATURE_GRANULARITY << "_out";
  S.matlab_input(filename.str().c_str());
  v.resize(S.row_dimension());
  for(unsigned int i=0; i< S.row_dimension(); i++)
   v[i] = S.get_entry( 0, i);
}


/* generates a time-stamp */
std::ostream& current_time(std::ostream& s)
{
  time_t rawtime = time(NULL);
  struct tm* timeinfo = localtime(&rawtime);
  char time_string[10];
  sprintf(time_string, "%02d:%02d:%02d", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
  s << time_string << "  ";
  return s;
}


/* MAIN */
int main()
{
  cout << "Testing simplified steepest descent algorithm in 2D with the biharmonic equation ..." << endl;

  const int DIM = 2;
  const int jmax = 4;


  //##############################  
  Matrix<double> A1(DIM,DIM);
  A1(0,0) = 2.;
  A1(1,1) = 1.0;
  Point<2> b1;
  b1[0] = -1.;
  b1[1] = -1.;
  AffineLinearMapping<2> affineP1(A1,b1);

  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 1.;
  A2(1,1) = 2.;
  Point<2> b2;
  b2[0] = -1.;
  b2[1] = -1.;
  AffineLinearMapping<2> affineP2(A2,b2);

  //##############################
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP1;
  charts[1] = &affineP2;
 
  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;
  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 2;
  bound_1[1] = 2;
  bound_1[2] = 2;
  bound_1[3] = 2;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 2;
  bound_2[1] = 2;
  bound_2[2] = 2;
  bound_2[3] = 2;

  bc[1] = bound_2;

  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bcT(2);

  //dual boundary conditions for first patch
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 0;
  bound_3[1] = 0;
  bound_3[2] = 0;
  bound_3[3] = 0;

  bcT[0] = bound_3;

  //dual boundary conditions for second patch
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 0;
  bound_4[1] = 0;
  bound_4[2] = 0;
  bound_4[3] = 0;
 
  bcT[1] = bound_4;

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  //finally a frame can be constructed
  //AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT, 5);
  cout << "Building frame from " << BASIS_NAME << " basis ..." << endl;
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, jmax);


  // setup biharmonic equation
  Vector<double> value(1);
  value[0] = 384;
  ConstantFunction<DIM> const_fun(value);

  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;

  CornerSingularityBiharmonic sing2D(origin, 0.5, 1.5);
  CornerSingularityBiharmonicRHS singRhs(origin, 0.5, 1.5);
  Functional<Basis1D, DIM> rhs(&singRhs, &frame);
  
  SimpleBiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&rhs, &frame,jmax, TrivialAffine);

  CachedProblem<SimpleBiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 7, 1.0/0.01);


  /* do iterative scheme */
  #if 1
  const double epsilon = 1.0e-6;

  InfiniteVector<double, Index> u_epsilon;

  current_time(cout);
  cout << "Performing iterative scheme ..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();
  cout << "here" <<endl;

  simplified_steepest_descent_SOLVE(problem, epsilon, u_epsilon, jmax);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  current_time(cout);
  cout << "Steepest descent done. Evaluating solution ..." << endl;

  u_epsilon.scale(&problem,-1);

  EvaluateFrame<Basis1D,2,2> evalObj;
  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, u_epsilon, true, PLOT_RESOLUTION);//expand in primal basis
  
  ostringstream filename_apprsol;
  filename_apprsol << "biharmonic_2D_steep_" << BASIS_NAME << "_jmax" << jmax << "_approx_sol.m";
  std::ofstream ofs51(filename_apprsol.str().c_str());
  matlab_output(ofs51,U);
  ofs51.close();

  current_time(cout);
  cout << "done plotting approximate solution, written file " << filename_apprsol.str() << endl;

  current_time(cout);
  cout << "Plotting pointwise error ..." << endl;

  Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, u_epsilon, sing2D, PLOT_RESOLUTION);
  //cout << "l_infty norm of error is " << linfty_norm(Error) << endl;

  ostringstream filename_error;
  filename_error << "biharmonic_2D_steep_" << BASIS_NAME << "_jmax" << jmax << "_error.m";
  std::ofstream ofs61(filename_error.str().c_str());
  matlab_output(ofs61,Error);
  ofs61.close();

  current_time(cout);
  cout << "done plotting pointwise error, written file " << filename_error.str() << endl;
#endif
  return 0;
}

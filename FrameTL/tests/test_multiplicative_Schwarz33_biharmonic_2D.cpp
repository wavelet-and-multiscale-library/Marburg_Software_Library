#define H_2_semi_norm_singularity 71.551 // preprocessed with matlab

#define BIHARMONIC

#define JMAX 8
#define TWO_D
#define SPARSE

#define OVERLAP 1.0

#define PATCHES 2

#define PRIMALORDER 3
#define DUALORDER   3

#define PRECOMP_RHS

#define PLOT_RESOLUTION 6

#define RHS_QUADRATURE_GRANULARITY 1

#define BASIS_NAME "P"

#include <fstream>
#include <iostream>
#include <time.h>

#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>

#include <interval/p_basis.h>
#include <interval/p_support.h>
#include <interval/p_evaluate.h>

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
#include <adaptive_multiplicative_Schwarz.h>

#include <map>

#include <cube/cube_indexplot.h>

typedef PBasis<PRIMALORDER,DUALORDER> Basis1D;
typedef AggregatedFrame<Basis1D,2,2> Frame2D;
typedef CubeBasis<Basis1D> Basis;
typedef MappedCubeBasis<Basis1D,2,2> MappedBasis;
typedef Frame2D::Index Index;
typedef CubeIndex<Basis1D,2,MappedCubeBasis<Basis1D,2,2> > CIndex;


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

  const int d = PRIMALORDER, dT = DUALORDER;
  const int DIM = 2;
  const int jmax = JMAX;


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
  
  SimpleBiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&rhs, &frame, jmax);

  //CachedProblem<SimpleBiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 7, 1.0/0.01);
  CachedProblemLocal<SimpleBiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 7, 1.0/0.01);

//   // plotting exact solution and right hand side
//   // grid resolution should be a multiple of three to get the re-entrant corners into the grid
//   Grid<1> grid1D(-1, 1, 512);
//   Grid<2> grid2D (grid1D, grid1D);
//   SampledMapping<2> mapping_exact (grid2D, sing2D);
//   SampledMapping<2> mapping_rhs (grid2D, singRhs);
//   std::ofstream of_exactsol("exact_solution_biharmL.m");
//   mapping_exact.matlab_output(of_exactsol);
//   of_exactsol.close();
//   std::ofstream of_rhs("rhs_biharmL.m");
//   mapping_rhs.matlab_output(of_rhs);
//   of_rhs.close();
//   // end plotting exact solution

//   abort();

  /* do iterative scheme */
#if 1
  const double epsilon = 1.0e-6;

  current_time(cout);
  cout << "Performing iterative scheme ..." << endl;
  clock_t tstart, tend;
  double time;
  tstart = clock();
  cout << "here" <<endl;


  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);
  MultSchw(problem, epsilon, approximations);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  current_time(cout);
  cout << "MultSchw done. Evaluating solution ..." << endl;

  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&problem,-1);
  
  EvaluateFrame<Basis1D,2,2> evalObj;
  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, approximations[frame.n_p()], true, PLOT_RESOLUTION);//expand in primal basis
  
  ostringstream filename_apprsol;
  filename_apprsol << "./ms_results2D_33_biharmL/biharmonic_2D_ms_" << BASIS_NAME << "_jmax" << jmax << "_approx_sol.m";
  std::ofstream ofs51(filename_apprsol.str().c_str());
  matlab_output(ofs51,U);
  ofs51.close();

  current_time(cout);
  cout << "done plotting approximate solution, written file " << filename_apprsol.str() << endl;

  current_time(cout);
  cout << "Plotting pointwise error ..." << endl;

  Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, approximations[frame.n_p()], sing2D, PLOT_RESOLUTION);
  //cout << "l_infty norm of error is " << linfty_norm(Error) << endl;

  ostringstream filename_error;
  filename_error << "./ms_results2D_33_biharmL/biharmonic_2D_ms_" << BASIS_NAME << "_jmax" << jmax << "_error.m";
  std::ofstream ofs61(filename_error.str().c_str());
  matlab_output(ofs61,Error);
  ofs61.close();

  current_time(cout);
  cout << "done plotting pointwise error, written file " << filename_error.str() << endl;

  for (int i = 0; i < frame.n_p(); i++) {
    cout << "plotting local approximation on patch " << i << endl;

    char filename3[128];
    sprintf(filename3, "%s%d%s%d%s%d%s", "./ms_results2D_33_biharmL/approx2D_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");

    U = evalObj.evaluate(frame, approximations[i], true, 6);//expand in primal basis
    std::ofstream ofsloc(filename3);
    matlab_output(ofsloc,U);
    //gnuplot_output(ofsloc,U);
    ofsloc.close();
  
  }

  cout << "potting sets of active wavelet indices..." << endl;
  Array1D<InfiniteVector<double, CIndex> > approximations_cube(frame.n_p());
  
  // convert indices to CubeIndices
  for (int i = 0; i < frame.n_p(); i++) {
    MappedBasis* mapped_basis = frame.bases()[i];
    std::ofstream plotstream;
    char filename4[128];
    sprintf(filename4, "%s%d%s%d%s%d%s", "./ms_results2D_33_biharmL/coefficient_plot_2D_patch_" , i , "_d" , d ,  "_dT", dT, ".m");
    plotstream.open(filename4);
    for (InfiniteVector<double, Index>::const_iterator it = approximations[frame.n_p()].begin(),
	   itend = approximations[frame.n_p()].end(); it != itend; ++it)
      if (it.index().p() == i) {
	approximations_cube[i].set_coefficient(CIndex(it.index().j(),it.index().e(),it.index().k(), mapped_basis),*it);
      }

    plot_indices_cube(frame.bases()[i], approximations_cube[i], jmax, plotstream, "jet", false, true);
    plotstream.close();

  }


#endif
  return 0;
}

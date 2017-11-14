#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0
#define _WAVELETTL_CDD1_VERBOSITY 0

#define H_1_semi_norm_singularity 1.6544 // preprocessed with high granularity

#define OVERLAP 1.

#define JMAX 7

#define RINGDOMAIN

#define PRIMALORDER 3
#define DUALORDER   3

#define PRECOMP_RHS
//#define PRECOMP_DIAG

//#define COMPUTECONSTANTS

#define SPARSE
//#define FULL
#define TWO_D

//#define PARALLEL

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#ifdef PARALLEL
#include <simple_elliptic_equation_parallel.h>
#else
#include <simple_elliptic_equation.h>
#endif
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <numerics/eigenvalues.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin/galerkin_utils.h>
#include <numerics/corner_singularity.h>
#include <frame_support.h>
#include <frame_index.h>
#ifdef PARALLEL
#include <adaptive_additive_Schwarz_parallel.h>
#else
#include <adaptive_additive_Schwarz.h>
#endif
#include <error_H_scale.h>
#include <galerkin/cached_problem.h>
#include <cube/cube_indexplot.h>
#include <poisson_2d_ring_testcase.h>

#ifdef PARALLEL
#include <parallel.h>
#endif


using std::cout;
using std::endl;

using FrameTL::FrameIndex;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
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

int main(int argc, char* argv[])
{
#ifdef PARALLEL
  MPI::Init (argc, argv);
  int rank, size;
  
  size = MPI::COMM_WORLD.Get_size ( );
  rank = MPI::COMM_WORLD.Get_rank ( );

  setup_coefficient_datatype();
#endif

  cout << "testing additive schwarz algorithm in 2D..." << endl;

  const int DIM = 2;

  const int jmax = JMAX;
  
  const int d = PRIMALORDER, dT = DUALORDER;

  //typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef CubeBasis<Basis1D> Basis;
  typedef MappedCubeBasis<Basis1D,2,2> MappedBasis;
  typedef Frame2D::Index Index;
  typedef CubeIndex<Basis1D,2,MappedCubeBasis<Basis1D,2,2> > CIndex;


  Matrix<double> A(DIM,DIM);
  A(0,0) = 3.0;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = -1.0;
  b[1] = -1.0;
  AffineLinearMapping<2> affineP(A,b);

  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 1.0;
  A2(1,1) = 3.0;
  Point<2> b2;
  b2[0] =  1.0;
  b2[1] = -1.0;
  AffineLinearMapping<2> affineP2(A2,b2);
 
  Matrix<double> A3(DIM,DIM);
  A3(0,0) = 3.0;
  A3(1,1) = 1.0;
  Point<2> b3;
  b3[0] = -1.0;
  b3[1] =  1.0;
  AffineLinearMapping<2> affineP3(A3,b3);

  Matrix<double> A4(DIM,DIM);
  A4(0,0) = 1.0;
  A4(1,1) = 3.0;
  Point<2> b4;
  b4[0] = -1.0;
  b4[1] = -1.0;
  AffineLinearMapping<2> affineP4(A4,b4);


  //##############################
  Array1D<Chart<DIM,DIM>* > charts(4);
  charts[0] = &affineP;
  charts[1] = &affineP2;
  charts[2] = &affineP3;
  charts[3] = &affineP4;


  SymmetricMatrix<bool> adj(4);

  // patch 0
  adj(0,0) = 1;
  adj(0,1) = 1;
  adj(0,3) = 1;

  // patch 1
  adj(1,0) = 1;
  adj(1,1) = 1;
  adj(1,2) = 1;

  // patch 2
  adj(2,1) = 1;
  adj(2,2) = 1;
  adj(2,3) = 1;

  // patch 3
  adj(3,0) = 1;
  adj(3,2) = 1;
  adj(3,3) = 1;

  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(4);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;//2
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 1;//2;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;//2;
  bound_2[2] = 1;
  bound_2[3] = 1;

  bc[1] = bound_2;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_3;
  bound_3[0] = 1;
  bound_3[1] = 1;//2;
  bound_3[2] = 1;
  bound_3[3] = 1;

  bc[2] = bound_3;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_4;
  bound_4[0] = 1;
  bound_4[1] = 1;//2;
  bound_4[2] = 1;
  bound_4[3] = 1;

  bc[3] = bound_4;



  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  //finally a frame can be constructed
  //AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT, 6);

  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, jmax);

  Vector<double> value(1);
  value[0] = 1;
  
  ConstantFunction<DIM> const_fun(value);
  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;

  Poisson_Solution_Ring sing2D;
  Poisson_RHS_Ring singRhs;
  Poisson_SolutionGradient_Ring singGrad();
  
//   // plotting exact solution and right hand side
//   // grid resolution should be a multiple of three to get the re-entrant corners into the grid
//   Grid<1> grid1D(-1, 2, 256);
//   Grid<2> grid2D (grid1D, grid1D);
//   SampledMapping<2> mapping_exact (grid2D, sing2D);
//   SampledMapping<2> mapping_rhs (grid2D, singRhs);
//   std::ofstream of_exactsol("exact_solution_ring.m");
//   mapping_exact.matlab_output(of_exactsol);
//   of_exactsol.close();
//   std::ofstream of_rhs("rhs_ring.m");
//   mapping_rhs.matlab_output(of_rhs);
//   of_rhs.close();
//   // end plotting exact solution

  
  PoissonBVP<DIM> poisson(&singRhs);
  //PoissonBVP<DIM> poisson(&const_fun);
  //PoissonBVP<DIM> poisson(&simple_sol_rhs);

  clock_t tstart, tend;
  double time;
#ifdef PARALLEL
  if (rank==0)
#endif
    tstart = clock();



  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);
  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);

  CachedProblemLocal<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
  discrete_poisson.set_norm_A(1.0);
  discrete_poisson.set_Ainv(1.0);

  const double epsilon = 1.0e-2;

  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);


  // ##########################################################################################
  // estimate extremal eigenvalues of local stiffness matrices and largest eigenvalue
  // of whole stiffness matrix
  
#ifdef COMPUTECONSTANTS
  set<Index> Lambda_0;
  set<Index> Lambda_1;
  set<Index> Lambda;
  for (FrameIndex<Basis1D,DIM,DIM> lambda = FrameTL::first_generator<Basis1D,DIM,DIM,Frame2D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,DIM,DIM,Frame2D>(&frame, jmax); ++lambda) {
    Lambda.insert(lambda);
    if (lambda.p() == 0)
      Lambda_0.insert(lambda);
    else {
      Lambda_1.insert(lambda);
    }
  }
  
  SparseMatrix<double> stiff;
  
  // starting vector for Power and Inverse Power Iteration
  Vector<double> x(Lambda_0.size()); x = 1;
  // number of iterations in Power and Inverse Power Iteration
  unsigned int iter= 0;  

  WaveletTL::setup_stiffness_matrix(problem, Lambda_0, stiff);

  cout << "computing smallest eigenvalue of stiffness matrix on patch 0" << endl;
  double lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "smallest eigenvalue of stiffness matrix on patch 0 is " << lmin << endl;

  cout << "computing largest eigenvalue of stiffness matrix on patch 0" << endl;
  double lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of stiffness matrix on patch 0 is " << lmax << endl;

  WaveletTL::setup_stiffness_matrix(problem, Lambda_1, stiff);

  x.resize(Lambda_1.size()); x = 1;
  cout << "computing smallest eigenvalue of stiffness matrix on patch 1" << endl;
  lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "smallest eigenvalue of stiffness matrix on patch 1 is " << lmin << endl;

  cout << "computing largest eigenvalue of stiffness matrix on patch 1" << endl;
  lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of stiffness matrix on patch 1 is " << lmax << endl;

  x.resize(Lambda.size()); x = 1;
  WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
  cout << "computing largest eigenvalue of whole stiffness matrix" << endl;
  lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  cout << "largest eigenvalue of whole stiffness matrix is " << lmax << endl;

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
  // ##########################################################################################

#ifdef PARLLEL
  if (rank==0) 
#endif
    cout << "adaptive algorithm started..." << endl;

  AddSchw(problem, epsilon, approximations);

#ifdef PARALLEL
  if (rank==0) {
#endif
    tend = clock();
    time = (double)(tend-tstart)/CLOCKS_PER_SEC;
    cout << "  ... done, time needed: " << time << " seconds" << endl;

    cout << "algorithm done" << endl;
  
    cout << "adaptive additive Schwarz done, generating output..." << endl;
 
    for (int i = 0; i <= frame.n_p(); i++)
      approximations[i].scale(&discrete_poisson,-1);
  
    EvaluateFrame<Basis1D,2,2> evalObj;

    Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, approximations[frame.n_p()], true, 6);//expand in primal basis
    cout << "...finished plotting global approximate solution" << endl;
    Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, approximations[frame.n_p()], sing2D, 6);
    //Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, approximations[frame.n_p()], simple_sol, 5);
    cout << "...finished plotting global error" << endl;

    char filename1[128];
#ifdef PARALLEL
    sprintf(filename1, "%s%d%s%d%s", "./as_results2D_33_par/approx2D_global_d", d, "_dT", dT, ".m");
#else
    sprintf(filename1, "%s%d%s%d%s", "./as_results2D_33/approx2D_global_d", d, "_dT", dT, ".m");
#endif
    std::ofstream ofs(filename1);
    matlab_output(ofs,U);
    //gnuplot_output(ofs,U);
    ofs.close();

    char filename2[128];
#ifdef PARALLEL
    sprintf(filename2, "%s%d%s%d%s", "./as_results2D_33_par/error2D_global_d", d, "_dT", dT, ".m");
#else
    sprintf(filename2, "%s%d%s%d%s", "./as_results2D_33/error2D_global_d", d, "_dT", dT, ".m");
#endif
    std::ofstream ofs1(filename2);
    matlab_output(ofs1,Error);
    //gnuplot_output(ofs1,Error);
    ofs1.close();


    for (int i = 0; i < frame.n_p(); i++) {
      cout << "plotting local approximation on patch " << i << endl;

      char filename3[128];
#ifdef PARALLEL
      sprintf(filename3, "%s%d%s%d%s%d%s", "./as_results2D_33_par/approx2D_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");
#else
      sprintf(filename3, "%s%d%s%d%s%d%s", "./as_results2D_33/approx2D_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");
#endif

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
#ifdef PARALLEL
      sprintf(filename4, "%s%d%s%d%s%d%s", "./as_results2D_33_par/coefficient_plot_2D_patch_" , i , "_d" , d ,  "_dT", dT, ".m");
#else
      sprintf(filename4, "%s%d%s%d%s%d%s", "./as_results2D_33/coefficient_plot_2D_patch_" , i , "_d" , d ,  "_dT", dT, ".m"); 
#endif
      plotstream.open(filename4);
      for (InfiniteVector<double, Index>::const_iterator it = approximations[frame.n_p()].begin(),
	     itend = approximations[frame.n_p()].end(); it != itend; ++it)
	if (it.index().p() == i) {
	  approximations_cube[i].set_coefficient(CIndex(it.index().j(),it.index().e(),it.index().k(), mapped_basis),*it);
	}

      plot_indices_cube(frame.bases()[i], approximations_cube[i], jmax, plotstream, "jet", false, true);
      plotstream.close();

    }
#ifdef PARALLEL
  }
#endif

  
#ifdef PARALLEL
  MPI::Finalize();
#endif

  return 0;

}


#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0

#define JMAX 8
#define TWO_D

//#define PRECOMP_RHS

#define PRIMALORDER 3
#define DUALORDER   3


#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <interval/spline_basis.h>
#include <simple_elliptic_equation.h>
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
#include <steepest_descent.h>
#include <galerkin/cached_problem.h>
#include <cube/cube_indexplot.h>

using std::cout;
using std::endl;

using FrameTL::FrameIndex;
using FrameTL::SimpleEllipticEquation;
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


int main()
{
  
  cout << "testing steepest descent algorithm in 2D..." << endl;

  const int DIM = 2;
  const int jmax = JMAX;

  const int d  = PRIMALORDER;
  const int dT = DUALORDER;

  //typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  
 
  //typedef SplineBasis<d,dT,P_construction> Basis1D;

  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  //typedef CubeBasis<Basis1D> Basis;
  typedef Frame2D::Index Index;
  typedef MappedCubeBasis<Basis1D,2,2> MappedBasis;
  typedef CubeIndex<Basis1D,2,MappedCubeBasis<Basis1D,2,2> > CIndex;


  EvaluateFrame<Basis1D,2,2> evalObj;

  //##############################  
  Matrix<double> A(DIM,DIM);
  // A(0,0) = 2.0;
  A(0,0) = 1.25;
  A(1,1) = 1.0;
  Point<2> b;
  //b[0] = -1.;
  b[0] = -0.25;
  b[1] = -1.;
  AffineLinearMapping<2> affineP(A,b);

  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 1.0;
  A2(1,1) = 2.0;
  Point<2> b2;
  b2[0] = -1.;
  b2[1] = -1.;
  AffineLinearMapping<2> affineP2(A2,b2);
  //##############################


  //##############################
  Array1D<Chart<DIM,DIM>* > charts(2);
  //charts[0] = &bezierP;
  charts[0] = &affineP;
  //charts[1] = &bezierP2;
  charts[1] = &affineP2;

  //charts[0] = &simpleaffine1;
  //charts[1] = &simpleaffine2;
 
  SymmetricMatrix<bool> adj(2);
  adj(0,0) = 1;
  adj(1,1) = 1;
  adj(1,0) = 1;
  adj(0,1) = 1;
  
  //to specify primal boundary the conditions
  Array1D<FixedArray1D<int,2*DIM> > bc(2);

  //primal boundary conditions for first patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_1;
  bound_1[0] = 1;
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 1; // d-1

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1; // d-1
  bound_2[2] = 1;
  bound_2[3] = 1;

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
  //AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT, jmax);
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, jmax);

  Vector<double> value(1);
  value[0] = 1;
  
  ConstantFunction<DIM> const_fun(value);

  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;
  CornerSingularity sing2D(origin, 0.5, 1.5);
  CornerSingularityRHS singRhs(origin, 0.5, 1.5);
  
  PoissonBVP<DIM> poisson(&singRhs);
  //PoissonBVP<DIM> poisson(&const_fun);

  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);

  //  L-shaped: (-1,1)x(-1,0) \cup (-1,0)x(-1,1), DSBasis
  
//     // (d,dT) = (2,2)
//     CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.0048, 1.0/0.146);
//     discrete_poisson.set_norm_A(5.0048);
//     // optimistic guess:
//     discrete_poisson.set_Ainv(1.0/0.146);

//   // (d,dT) = (3,5)
//   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 11.7375, 1.0/0.146);
//   discrete_poisson.set_norm_A(11.7375);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);

//   // (d,dT) = (3,5) all boundary conditions of order 1
//   CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.0022, 1.0/0.146);
//   discrete_poisson.set_norm_A(8.0022);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);


//   // (d,dT) = (3,3)
//   //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 7.6336, 1.0/0.146);
//   CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 7.6336, 1.0/0.146);
//   discrete_poisson.set_norm_A(7.6336);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);



//   // (d,dT) = (4,6)
//   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 19.1803, 1.0/0.146);
//   discrete_poisson.set_norm_A(19.1803);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);


  //  L-shaped: (-1,1)x(-1,0) \cup (-1,0)x(-1,1), PBasis

//   // (d,dT) = (2,2)
//   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.0225, 1.0/0.146);
//   discrete_poisson.set_norm_A(5.0225);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);

//   // (d,dT) = (3,5)
//   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
//   discrete_poisson.set_norm_A(8.3898);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);
  
  // (d,dT) = (3,3)
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
  discrete_poisson.set_norm_A(1.0);
  // optimistic guess:
  discrete_poisson.set_Ainv(1.0);
//   CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson,7.2346, 1.0/0.146);
//   discrete_poisson.set_norm_A(7.2346);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);


//   // (d,dT) = (4,6)
//   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 14.6912, 1.0/0.146);
//   discrete_poisson.set_norm_A(14.6912);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);

  //   discrete_poisson.set_norm_A(6.99614824235842);
  //   discrete_poisson.set_Ainv(1.0/0.1);
  // CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.0048, 1.0/0.01);
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 6.99614824235842, 1.0/0.1);


  const double epsilon = 1.0e-6;

  //InfiniteVector<double, Index> u_epsilon;
  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);

  clock_t tstart, tend;
  double time;
  tstart = clock();

  //steepest_descent_SOLVE(problem, epsilon, u_epsilon, approximations);
  steepest_descent_SOLVE(problem, epsilon, approximations);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "steepest descent done" << endl;

  //discrete_poisson.rescale(u_epsilon,-1);
  //u_epsilon.scale(&problem,-1);
  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&discrete_poisson,-1);

#if 1

  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, approximations[frame.n_p()], true, 6);//expand in primal basis
  cout << "done plotting approximate solution" << endl;

  Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, approximations[frame.n_p()], sing2D, 6);
  cout << "done plotting pointwise error" << endl;

  std::ofstream ofs5("./sd_results2D_33/approx_sol_steep_2D_out.m");
  matlab_output(ofs5,U);
  ofs5.close();

  std::ofstream ofs6("./sd_results2D_33/error_steep_2D_out.m");
  matlab_output(ofs6,Error);
  ofs6.close();

  for (int i = 0; i < frame.n_p(); i++) {
    cout << "plotting local approximation on patch " << i << endl;

    char filename3[128];
    sprintf(filename3, "%s%d%s%d%s%d%s", "./sd_results2D_33/approx2Dsteep_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");

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
    char filename4[50];
    sprintf(filename4, "%s%d%s%d%s%d%s", "./sd_results2D_33/coefficient_plot_2D_patch_" , i , "_d" , d ,  "_dT", dT, ".m");
    plotstream.open(filename4);
    for (InfiniteVector<double, Index>::const_iterator it = approximations[frame.n_p()].begin(),
	   itend = approximations[frame.n_p()].end(); it != itend; ++it)
      if (it.index().p() == i) {
	approximations_cube[i].set_coefficient(CIndex(it.index().j(),it.index().e(),it.index().k(), mapped_basis),*it);
      }

    plot_indices(frame.bases()[i], approximations_cube[i], jmax, plotstream, "jet", false, true);
    plotstream.close();

  }

#endif
  //  problem.add_level(frame.first_generator(3),u_epsilon,3,1.);

  return 0;

}

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0

#define JMAX 6
#define TWO_D


#define PRIMALORDER 4
#define DUALORDER   6

#define CONSTFUN
#define PROJ

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
#include <richardson.h>
#include <richardson_projector_galerkin.h>
#include <galerkin/cached_problem.h>
#include <utils/plot_tools.h>
#include <interval/i_indexplot.h>

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

//using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
  
  cout << "testing Richardson iteration with projection step in 2D..." << endl;

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

  EvaluateFrame<Basis1D,2,2> evalObj;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 2.0;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = -1.;
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
  bound_1[3] = 1;//2

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;//2
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
  value[0] = 1.0;
  
  ConstantFunction<DIM> const_fun(value);

  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;
  CornerSingularity sing2D(origin, 0.5, 1.5);
//  CornerSingularityRHS singRhs(origin, 0.5, 1.5);
  
  // PoissonBVP<DIM> poisson(&singRhs);
  PoissonBVP<DIM> poisson(&const_fun);

  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);

  //  L-shaped: (-1,1)x(-1,0) \cup (-1,0)x(-1,1), DSBasis

  // (d,dT) = (3,3)
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
  discrete_poisson.set_norm_A(1.0);
  // optimistic guess:
  discrete_poisson.set_Ainv(1.0);

  const double epsilon = 1e-5;

  InfiniteVector<double, Index> u_epsilon;
  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);

  clock_t tstart, tend;
  double time;
  tstart = clock();

#ifdef PROJ
  richardson_SOLVE_projector_galerkin(problem, epsilon, u_epsilon, approximations);
#else
  richardson_SOLVE(problem, epsilon, u_epsilon, approximations);
#endif

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "Richardson done" << endl;

  //discrete_poisson.rescale(u_epsilon,-1);
  u_epsilon.scale(&problem,-1);
#if 1
  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, u_epsilon, true, 6);//expand in primal basis
  
  cout << "done plotting approximate solution" << endl;

  Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, u_epsilon, sing2D, 6);

  cout << "done plotting pointwise error" << endl;

  std::ofstream ofs5("./Richardson_Proj_2D/approx_sol_rich_2D_out.m");
  matlab_output(ofs5,U);
  ofs5.close();

  std::ofstream ofs6("./Richardson_Proj_2D/error_rich_2D_out.m");
  matlab_output(ofs6,Error);
  ofs6.close();

  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&discrete_poisson,-1);


  for (int i = 0; i < frame.n_p(); i++) {
    cout << "plotting local approximation on patch " << i << endl;

    char filename3[150];
    sprintf(filename3, "%s%d%s%d%s%d%s", "./Richardson_Proj_2D/approx2Drich_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");

    U = evalObj.evaluate(frame, approximations[i], true, 6);//expand in primal basis
    std::ofstream ofsloc(filename3);
    matlab_output(ofsloc,U);
    //gnuplot_output(ofsloc,U);
    ofsloc.close();

  }

  typedef Basis1D::Index Index1D;

  FixedArray1D<InfiniteVector<double, Index1D>, 2> indices;

  InfiniteVector<double, Index>::const_iterator it = approximations[frame.n_p()].begin();
  for (; it!= approximations[frame.n_p()].end(); ++it) {
    //cout << *it << endl;
    Index ind(it.index());
    //cout << "level = " << ind.j() << endl;
    indices[ind.p()].set_coefficient(Index1D(ind.j(),ind.e()[0],ind.k()[0],
					     frame.bases()[0]->bases()[0]), *it);

    //cout << log10(fabs(*it)) << endl;
  }

  std::ofstream ofs7("./Richardson_Proj_2D/indices_patch_0.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[0]->bases()[0], indices[0], JMAX, ofs7, "jet", true, -16);

  std::ofstream ofs8("./Richardson_Proj_2D/indices_patch_1.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[1]->bases()[0], indices[1], JMAX, ofs8, "jet", true, -16);
  // compute infinite vectors of 1D indices, one for each patch
  // and plot them


#endif
  //  problem.add_level(frame.first_generator(3),u_epsilon,3,1.);

  return 0;

}

#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0

#define JMAX 15
#define ONE_D

// define macro LOCALGALERKIN if Galerkin-method shall be used for local solvers
#define LOCALGALERKIN

// if LOCALGALERKIN is defined: define macro PROJECTOR if projection strategy shall be used
#define PROJECTOR

#define WAVELETTL_CDD1_VERBOSITY 0

#define PRIMALORDER 3
#define DUALORDER   3

#include <fstream> 
#include <iostream>
#include <time.h>
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <elliptic_equation.h>
#include <simple_elliptic_equation.h>
#include <biharmonic_equation.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin/galerkin_utils.h>
#include <frame_support.h>
#include <frame_index.h>
#include <richardson.h>
#include <richardson_projector.h>
#include <richardson_projector_galerkin.h>
#include <galerkin/cached_problem.h>
#include <utils/plot_tools.h>
#include <interval/i_indexplot.h>

using std::cout;
using std::endl;

//using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;


int main()
{
 
  cout << "testing Richardson iteration with projection step in 1D..." << endl;
  
  const int jmax = JMAX;
  const int DIM = 1;

  const int d = PRIMALORDER, dT = DUALORDER;

    typedef PBasis<d,dT> Basis1D;
  typedef AggregatedFrame<Basis1D,1,1> Frame1D;
  typedef CubeBasis<Basis1D,1> IntervalBasis;
  typedef Frame1D::Index Index;


  //##############################
  Matrix<double> A(DIM,DIM);
  A(0,0) = 0.7;
  Point<1> b;
  b[0] = 0.;
  AffineLinearMapping<1> affineP(A,b);
  
  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 0.7;
  Point<1> b2;
  b2[0] = 1-A2.get_entry(0,0);
  AffineLinearMapping<1> affineP2(A2,b2);


  //##############################
  
  Array1D<Chart<DIM,DIM>* > charts(2);
  charts[0] = &affineP;
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
  bound_1[0] = 1;
  bound_1[1] = d-1;
  
  bc[0] = bound_1;
  
  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = d-1;
  bound_2[1] = 1;
  
  bc[1] = bound_2;
  
  //to specify primal boundary the conditions
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

  Atlas<DIM,DIM> interval_domain(charts,adj);  
  cout << interval_domain << endl;

  //finally a frame can be constructed
  //Frame1D frame(&Lshaped, bc, bcT, jmax);
  Frame1D frame(&interval_domain, bc, jmax);

  Singularity1D_RHS_2<double> sing1D;
  Singularity1D_2<double> exactSolution1D;
 
  PoissonBVP<DIM> poisson(&sing1D);


  // models the discrete problem
  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);

  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 1.0, 1.0);
  discrete_poisson.set_norm_A(1.0);
  // optimistic guess:
  discrete_poisson.set_Ainv(1.0);

  const double epsilon = 1.0e-4;
  InfiniteVector<double, Index> u_epsilon;
  u_epsilon.clear();


  clock_t tstart, tend;
  double time;
  tstart = clock();

  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);

  // call solver
#ifdef LOCALGALERKIN
  richardson_SOLVE_projector_galerkin(problem, epsilon, u_epsilon, approximations);
#else
  richardson_SOLVE_projector(problem, epsilon, u_epsilon, approximations);
#endif


  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "Richardson iteration done" << endl;

  u_epsilon.scale(&discrete_poisson,-1);
  cout << "u_epsilon.size(): " << u_epsilon.size() << endl;

  EvaluateFrame<Basis1D,1,1> evalObj;

  Array1D<SampledMapping<1> > U = evalObj.evaluate(frame, u_epsilon, true, 12);//expand in primal basis
  cout << "...finished plotting approximate solution" << endl;
  Array1D<SampledMapping<1> > Error = evalObj.evaluate_difference(frame, u_epsilon, exactSolution1D, 12);
  cout << "...finished plotting error" << endl;
  
  std::ofstream ofs5("./Richardson_Proj_1D/approx_sol_rich_1D_out.m");
  matlab_output(ofs5,U);
  ofs5.close();

  std::ofstream ofs6("./Richardson_Proj_1D/error_rich_1D_out.m");
  matlab_output(ofs6,Error);
  ofs6.close();

  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&discrete_poisson,-1);
  
  for (int i = 0; i < frame.n_p(); i++) {
    cout << "plotting local approximation on patch " << i << endl;

    char filename3[128];
    sprintf(filename3, "%s%d%s%d%s%d%s", "./Richardson_Proj_1D/approx1Drich_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");

    U = evalObj.evaluate(frame, approximations[i], true, 12);//expand in primal basis
    std::ofstream ofsloc(filename3);
    matlab_output(ofsloc,U);
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

  std::ofstream ofs7("./Richardson_Proj_1D/indices_patch_0.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[0]->bases()[0], indices[0], JMAX, ofs7, "jet", true, -16);

  std::ofstream ofs8("./Richardson_Proj_1D/indices_patch_1.m");
  WaveletTL::plot_indices<Basis1D>(frame.bases()[1]->bases()[0], indices[1], JMAX, ofs8, "jet", true, -16);
  // compute infinite vectors of 1D indices, one for each patch
  // and plot them

  return 0;

}

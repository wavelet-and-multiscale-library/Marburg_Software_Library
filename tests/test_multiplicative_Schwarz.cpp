#define _WAVELETTL_GALERKINUTILS_VERBOSITY 0

#define OVERLAP 1.

//#define SPARSE
#define FULL
#define TWO_D

#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <elliptic_equation.h>
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
#include <adaptive_multiplicative_Schwarz.h>
//#include <multiplicative_Schwarz.h>
//#include <additive_Schwarz.h>
//#include <additive_Schwarz_SD.h>
#include <galerkin/cached_problem.h>

using std::cout;
using std::endl;

using FrameTL::FrameIndex;
using FrameTL::EllipticEquation;
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
 
  cout << "testing multiplicative schwarz algorithm in 2D..." << endl;

  const int DIM = 2;

  const int jmax = 5;
  
  const int d = 3, dT = 3;

  //typedef DSBasis<d,dT> Basis1D;
  typedef PBasis<d,dT> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef CubeBasis<Basis1D> Basis;
  typedef Frame2D::Index Index;

//   //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 1. + OVERLAP;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = -OVERLAP;
  b[1] = -1.;
  AffineLinearMapping<2> affineP(A,b);

  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 1.0;
  A2(1,1) = 2.0;
  Point<2> b2;
  b2[0] = -1.0;
  b2[1] = -1.0;
  AffineLinearMapping<2> affineP2(A2,b2);


//   LinearBezierMapping bezierP(Point<2>(-1.,-1.),Point<2>(-1.,1.),
//  			      Point<2>(0.,-1.), Point<2>(0.,1.));
  
//   LinearBezierMapping bezierP2(Point<2>(-1.,-1.),Point<2>(-1.,0.),
// 			       Point<2>(1.,-1.), Point<2>(1.,0.));
 

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
  bound_1[0] = 1;//2
  bound_1[1] = 1;
  bound_1[2] = 1;
  bound_1[3] = 2;//2;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 2;//2;
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
  //AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, bcT, 6);
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


  clock_t tstart, tend;
  double time;
  tstart = clock();


  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, TrivialAffine);
  //EllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, Composite);
  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);

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

  // (d,dT) = (3,5)
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  //CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  CachedProblemLocal<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  discrete_poisson.set_norm_A(8.3898);
  // optimistic guess:
  discrete_poisson.set_Ainv(1.0/0.146);

//   // (d,dT) = (4,6)
//   CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 14.6912, 1.0/0.146);
//   discrete_poisson.set_norm_A(14.6912);
//   // optimistic guess:
//   discrete_poisson.set_Ainv(1.0/0.146);

  //   discrete_poisson.set_norm_A(6.99614824235842);
  //   discrete_poisson.set_Ainv(1.0/0.1);
  // CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 5.0048, 1.0/0.01);
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 6.99614824235842, 1.0/0.1);



  const double epsilon = 0.01;

  Array1D<InfiniteVector<double, Index> > approximations(frame.n_p()+1);

//   set<Index> Lambda;
//   for (FrameIndex<Basis1D,2,2> lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
//        lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, jmax); ++lambda) {
//     Lambda.insert(lambda);
//     //cout << lambda << endl;
//   }
  
//   SparseMatrix<double> stiff;
  
//   WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);


  adaptive_multiplicative_Schwarz_SOLVE(problem, epsilon, approximations);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "algorithm done" << endl;

  cout << "adaptive multiplicative Schwarz done, generating output..." << endl;
 
  for (int i = 0; i <= frame.n_p(); i++)
    approximations[i].scale(&discrete_poisson,-1);
  
  EvaluateFrame<Basis1D,2,2> evalObj;

  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, approximations[frame.n_p()], true, 5);//expand in primal basis
  cout << "...finished plotting global approximate solution" << endl;
  Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, approximations[frame.n_p()], sing2D, 5);
  cout << "...finished plotting global error" << endl;

  char filename1[50];
  sprintf(filename1, "%s%d%s%d%s", "approx2D_global_d", d, "_dT", dT, ".m");
  std::ofstream ofs(filename1);
  matlab_output(ofs,U);
  ofs.close();

  char filename2[50];
  sprintf(filename2, "%s%d%s%d%s", "error2D_global_d", d, "_dT", dT, ".m");
  std::ofstream ofs1(filename2);
  matlab_output(ofs1,Error);
  ofs1.close();


  for (int i = 0; i < frame.n_p(); i++) {
    cout << "plotting local approximation on patch " << i << endl;

    char filename3[50];
    sprintf(filename3, "%s%d%s%d%s%d%s", "approx2D_local_on_patch_" , i , "_d" , d ,  "_dT", dT, ".m");

    U = evalObj.evaluate(frame, approximations[i], true, 5);//expand in primal basis
    std::ofstream ofsloc(filename3);
    matlab_output(ofsloc,U);
    ofsloc.close();
  
  }
  return 0;

}


#include <fstream>
#include <iostream>
#include <time.h> 
#include <interval/ds_basis.h>
#include <interval/p_basis.h>
#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <numerics/eigenvalues.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
#include <galerkin/galerkin_utils.h>
#include <frame_support.h>
#include <frame_index.h>
#include <steepest_descent1.h>
#include <galerkin/cached_problem.h>
#include <map>
#include <numerics/corner_singularity_biharmonic.h>
#include <biharmonic_equation.h>



using std::cout;
using std::endl;

using FrameTL::FrameIndex;
using FrameTL::BiharmonicEquation;
using FrameTL::EvaluateFrame;
using FrameTL::AggregatedFrame;
using MathTL::BiharmonicBVP;
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


int main()
{
  
  cout << "testing steepest descent algorithm in 2D..." << endl;

 const int DIM = 2;
  const int jmax = 5;

  //typedef DSBasis<4,6> Basis1D;
  typedef PBasis<3,3> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  typedef CubeBasis<Basis1D> Basis;
  typedef Frame2D::Index Index;

  EvaluateFrame<Basis1D,2,2> evalObj;

  //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 2.;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = -1.;
  b[1] = -1.;
  AffineLinearMapping<2> affineP(A,b);

  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 1.;
  A2(1,1) = 2.;
  Point<2> b2;
  b2[0] = -1.;
  b2[1] = -1.;
  AffineLinearMapping<2> affineP2(A2,b2);
  //##############################

  //##############################
  LinearBezierMapping bezierP(Point<2>(-1.,-1.),Point<2>(-1.,1.),
 			      Point<2>(0.,-1.), Point<2>(0.,1.));
  
  LinearBezierMapping bezierP2(Point<2>(-1.,-1.),Point<2>(-1.,0.),
			       Point<2>(1.,-1.), Point<2>(1.,0.));
 


  FixedArray1D<double,2> A3;
  A3[0] = 1.;
  A3[1] = 2.;
  SimpleAffineLinearMapping<2> simpleaffine1(A3,b);
  
  FixedArray1D<double,2> A4;
  A4[0] = 2.;
  A4[1] = 1.;
  SimpleAffineLinearMapping<2> simpleaffine2(A4,b2);


  //##############################
  Array1D<Chart<DIM,DIM>* > charts(2);
  //charts[0] = &bezierP;
  charts[0] = &affineP;
  //charts[1] = &bezierP2;
  charts[1] = &affineP2;
 
  //   charts[0] = &simpleaffine1;
  //   charts[1] = &simpleaffine2;

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

  Atlas<DIM,DIM> Lshaped(charts,adj);  
  cout << Lshaped << endl;

  //finally a frame can be constructed
  AggregatedFrame<Basis1D, DIM, DIM> frame(&Lshaped, bc, jmax);

  Vector<double> value(1);
  value[0] = 1;
  
  ConstantFunction<DIM> const_fun(value);
  Point<2> origin;
  origin[0] = 0.0;
  origin[1] = 0.0;


  CornerSingularityBiharmonic sing2D(origin, 0.5, 1.5);
  CornerSingularityBiharmonicRHS singRhs(origin, 0.5, 1.5);

  BiharmonicBVP<DIM> biharmonic(&singRhs);
  
  BiharmonicEquation<Basis1D,DIM> discrete_biharmonic(&biharmonic, &frame,jmax, TrivialAffine);  

  CachedProblem<BiharmonicEquation<Basis1D,DIM> > problem(&discrete_biharmonic, 5.0048, 1.0/0.01);

  const double epsilon = 1.0e-6;

  InfiniteVector<double, Index> u_epsilon;

  clock_t tstart, tend;
  double time;
  tstart = clock();

  steepest_descent1_SOLVE(problem, epsilon, u_epsilon, jmax);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

  cout << "steepest descent done" << endl;

  u_epsilon.scale(&problem,-1);

  Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, u_epsilon, true, 6);//expand in primal basis
  
  cout << "done plotting approximate solution" << endl;

  Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, u_epsilon, sing2D, 6);

  cout << "done plotting pointwise error" << endl;

  std::ofstream ofs5("approx_sol_steep_2D_out.m");
  matlab_output(ofs5,U);
  ofs5.close();

  std::ofstream ofs6("error_steep_2D_out.m");
  matlab_output(ofs6,Error);
  ofs6.close();

  //  problem.add_level(frame.first_generator(3),u_epsilon,3,1.);

  return 0;

}

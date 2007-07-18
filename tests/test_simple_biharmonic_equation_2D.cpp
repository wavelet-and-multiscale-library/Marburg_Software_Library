#define MAX_LOOPS 10000

#define BASIS_S

//#define USED_SAVED_RHS
//#define SAVE_RHS

#define RHS_QUADRATURE_GRANULARITY 4

#include <fstream>
#include <iostream>
#include <time.h>

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

#include <numerics/corner_singularity_biharmonic.h>
#include <simple_biharmonic_equation.h>

#include <algebra/sparse_matrix.h>
#include <algebra/infinite_vector.h>
#include <numerics/iteratsolv.h>
#include <frame_evaluate.h>
#include <cube/cube_basis.h>
#include <cube/cube_index.h>
//#include <galerkin_frame_utils.h>
#include <galerkin/galerkin_utils.h>
#include <galerkin/cached_problem.h>
#include <frame_support.h>
#include <numerics/eigenvalues.h>
//#include <biharmonic_rhs.h>
//#include <numerics/bvp.h>

using std::cout;
using std::endl;

using FrameTL::EvaluateFrame;
//using FrameTL::EllipticEquation;
using FrameTL::SimpleBiharmonicEquation;
using FrameTL::AggregatedFrame;
using MathTL::ConstantFunction;
using MathTL::CornerSingularityBiharmonicRHS;
using MathTL::SparseMatrix;
using MathTL::InfiniteVector;
using WaveletTL::CubeBasis;
using WaveletTL::CubeIndex;

using namespace std;
using namespace FrameTL;
using namespace MathTL;
using namespace WaveletTL;

//#define TWO_DIMENSIONS

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


std::ostream& current_time(std::ostream& s)
{
  time_t rawtime = time(NULL);
  struct tm* timeinfo = localtime(&rawtime);
  char time_string[10];
  sprintf(time_string, "%02d:%02d:%02d", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
  s << time_string << "  ";
  return s;
}


int main()
{
  
  cout << "Testing class SimpleBiharmonicEquation in 2D ..." << endl;
  
  const int DIM = 2;
  const int jmax = 3;
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

//   double tmp = 0.0;
//   int c = 0;
//   int d = 0;
  
  cout.precision(12);

  //############### 2D galerkin scheme test ##################
  //#if 1

  set<Index> Lambda;
  for (FrameIndex<Basis1D,2,2> lambda = FrameTL::first_generator<Basis1D,2,2,Frame2D>(&frame, frame.j0());
       lambda <= FrameTL::last_wavelet<Basis1D,2,2,Frame2D>(&frame, jmax); ++lambda) {
    Lambda.insert(lambda);
    //cout << lambda << endl;
  }

  current_time(cout);
  cout << "Setting up full right hand side ..." << endl;
  Vector<double> rh;
  #ifdef USE_SAVED_RHS
  Rhs_input(rh, jmax);
  #else
  WaveletTL::setup_righthand_side(discrete_biharmonic, Lambda, rh);
   #ifdef SAVE_RHS
   cout << "Saving right hand side to file ..." << endl;
   Rhs_output(rh, jmax); 
   #endif
  #endif
  //cout << rh << endl;
  cout << "l_2-norm of RHS (with quadrature granularity " << RHS_QUADRATURE_GRANULARITY << "): " << sqrt(rh*rh) << endl;
#if 1
  current_time(cout);
  cout << "Setting up full stiffness matrix ..." << endl;
  SparseMatrix<double> stiff;

  clock_t tstart, tend;
  double time;
  tstart = clock();

  WaveletTL::setup_stiffness_matrix(discrete_biharmonic, Lambda, stiff);
  //WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);

  tend = clock();
  time = (double)(tend-tstart)/CLOCKS_PER_SEC;
  cout << "  ... done, time needed: " << time << " seconds" << endl;

#ifdef SAVE_STIFFNESS_MATRIX
  cout << "Saving stiffness matrix to file ..." << endl;
  ostringstream filename_stiff;
  filename_stiff << "biharmonic_stiff_2D_" << BASIS_NAME << "_out";
  stiff.matlab_output(filename_stiff.str().c_str(), "stiff",1);
#endif
  
//   unsigned int k;
//   double lambdamin, lambdamax;
//   LanczosIteration(stiff, 1e-6,
// 			lambdamin, lambdamax,
// 		   200, k);

//   cout << "lambdamax=" << lambdamax << "  lambdamin=" << lambdamin << endl;


  //unsigned int iter= 0;
  //ctor<double> x(Lambda.size()); x = 1;
  //double lmax = PowerIteration(stiff, x, 0.01, 1000, iter);
  double lmax = 1;
  //ut << "lmax = " << lmax << endl;

  //x = 1;
  //double lmin = InversePowerIteration(stiff, x, 0.01, 1000, iter);

  current_time(cout);
  cout << "Performing iterative scheme to solve projected problem ..." << endl;
  
  Vector<double> xk(Lambda.size()), err(Lambda.size());
  xk = 0;
  
  //CG(stiff, rh, xk, 1.0e-6, 100, iter);
  //cout << "CG iterations needed: "  << iter << endl;
  //Richardson(stiff, rh, xk, 2. / lmax - 0.01, 0.0001, 2000, iter);
  
  double alpha_n = 2. / lmax - 0.001;
  
  Vector<double> resid(xk.size());
  Vector<double> help(xk.size());

//   for (int i = 0; i < 1500;i++) {
//     stiff.apply(xk,help);
//     resid = rh - help;
//     cout << "i = " << i << " " << sqrt(resid*resid) << endl;
//     stiff.apply(resid,help);
//     alpha_n = (resid * resid) * (1.0 / (resid * help));
//     resid *= alpha_n;
//     xk = xk + resid;
//   }




//  cout << "performing output ..." << endl;
  
   InfiniteVector<double,Frame2D::Index> u;
   unsigned int i = 0;
//   for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
//     u.set_coefficient(*it, xk[i]);
  
//   u.scale(&discrete_biharmonic,-1);
 
//   EvaluateFrame<Basis1D,2,2> evalObj;
//   Array1D<SampledMapping<2> > U = evalObj.evaluate(frame, u, true, jmax);//expand in primal basis
  
//   Array1D<SampledMapping<2> > Error = evalObj.evaluate_difference(frame, u, sing2D, 6);

//   std::ofstream ofs6("error_2D_out_jmax6_steps1500_b4.m");
//   matlab_output(ofs6,Error);
//   ofs6.close();

//   std::ofstream ofs5("approx_solution_out_6_1500_b4.m");
//   matlab_output(ofs5,U);
//   ofs5.close();

  xk=0;
  alpha_n = 2. / lmax - 0.001;
  double min=1;
  unsigned int z = 0;
  double res_norm;
  for (unsigned int i = 0; i < MAX_LOOPS; i++) {
    stiff.apply(xk,help);
    resid = rh - help;
    res_norm = sqrt(resid*resid);
    current_time(cout);
    cout << "loop " << i << ": res_norm = " << res_norm << endl;
    if (res_norm < min) {
      min = res_norm;
      z = i;
    }
    stiff.apply(resid,help);
    alpha_n = (resid * resid) * (1.0 / (resid * help));
    resid *= alpha_n;
    xk = xk + resid;
  }


  current_time(cout);
  cout << "Performing output ..." << endl;

  i = 0;
  for (set<Index>::const_iterator it = Lambda.begin(); it != Lambda.end(); ++it, ++i)
    u.set_coefficient(*it, xk[i]);
  
  u.scale(&discrete_biharmonic,-1);
 
  EvaluateFrame<Basis1D,2,2> evalObj;
  Array1D<SampledMapping<2> > U1 = evalObj.evaluate(frame, u, true, jmax);//expand in primal basis
  
  Array1D<SampledMapping<2> > Error1 = evalObj.evaluate_difference(frame, u, sing2D, 6);

  ostringstream filename_error;
  filename_error << "biharmonic_2D_" << BASIS_NAME << "_jmax" << jmax << "_error.m";
  std::ofstream ofs61(filename_error.str().c_str());
  matlab_output(ofs61,Error1);
  ofs61.close();

  ostringstream filename_apprsol;
  filename_apprsol << "biharmonic_2D_" << BASIS_NAME << "_jmax" << jmax << "_approx_solution_out.m";
  std::ofstream ofs51(filename_apprsol.str().c_str());
  matlab_output(ofs51,U1);
  ofs51.close();
  cout << "Minimum:" << min << " im Schritt " << z << endl;
#endif
  

  return 0;
}

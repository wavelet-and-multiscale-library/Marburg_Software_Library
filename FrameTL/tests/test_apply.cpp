#define _WAVELETTL_GALERKINUTILS_VERBOSITY 1
#define PARALLEL 0
#define _DIM 2

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
#include <multiplicative_Schwarz.h>
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

  const int jmax = 10;

  //typedef DSBasis<4,6> Basis1D;
  typedef PBasis<2,2> Basis1D;
  typedef AggregatedFrame<Basis1D,2,2> Frame2D;
  //typedef CubeBasis<Basis1D> Basis;
  typedef Frame2D::Index Index;

  //EvaluateFrame<Basis1D,2,2> evalObj;

//   //##############################  
  Matrix<double> A(DIM,DIM);
  A(0,0) = 2.;
  A(1,1) = 1.0;
  Point<2> b;
  b[0] = -1.;
  b[1] = -1.;
  AffineLinearMapping<2> affineP(A,b);

  Matrix<double> A2(DIM,DIM);
  A2(0,0) = 1.0;
  A2(1,1) = 2.0;
  Point<2> b2;
  b2[0] = -1.0;
  b2[1] = -1.0;
  AffineLinearMapping<2> affineP2(A2,b2);

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
  bound_1[3] = 1;//2;

  bc[0] = bound_1;

  //primal boundary conditions for second patch: all Dirichlet
  FixedArray1D<int,2*DIM> bound_2;
  bound_2[0] = 1;
  bound_2[1] = 1;//2;
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
  
  //PoissonBVP<DIM> poisson(&singRhs);
  PoissonBVP<DIM> poisson(&const_fun);


  SimpleEllipticEquation<Basis1D,DIM> discrete_poisson(&poisson, &frame, jmax);

  // (d,dT) = (3,5)
  //CachedProblem<EllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  CachedProblem<SimpleEllipticEquation<Basis1D,DIM> > problem(&discrete_poisson, 8.3898, 1.0/0.146);
  discrete_poisson.set_norm_A(8.3898);
  // optimistic guess:
  discrete_poisson.set_Ainv(1.0/0.146);


  cout.precision(12);
  
  InfiniteVector<double,Frame2D::Index> rhs;
  Vector<double> v(225);
  
//   set<Index> Lambda;
//   for (FrameIndex<Basis1D,DIM,DIM> lambda = FrameTL::first_generator<Basis1D,DIM,DIM,Frame2D>(&frame, frame.j0());
//        lambda <= FrameTL::last_wavelet<Basis1D,DIM,DIM,Frame2D>(&frame, jmax); ++lambda) {
//     Lambda.insert(lambda);
//   }
  
//   cout << "setting up full right hand side..." << endl;
//   Vector<double> rh;
//   WaveletTL::setup_righthand_side(discrete_poisson, Lambda, rh);
//   //cout << rh << endl;

  cout << "setting up full stiffness matrix..." << endl;
  SparseMatrix<double> stiff;

  //WaveletTL::setup_stiffness_matrix(problem, Lambda, stiff);
  //cout << "NUMBER OF NONZERO ENTRIES = " << stiff.size() << endl;

  InfiniteVector<double, Index> w, f, w_res1, w_res2, w_eps, w_exact, error, sparse_f;
  
//   // copy rh to w
//   for (unsigned int i = 0; i < rh.size(); i++) {
//     if (rh[i] != 0.) {
//       Index ind = *problem.basis().get_wavelet(i);
//       f.set_coefficient(ind, rh[i]);
//     }
//   }
  
  cout << "number of nonzeros in right hand side " << w.size() << endl;


  
//   //APPLY_TEST(problem, f, 0.0001, w_res1, jmax, St04a);
//   APPLY(problem, f, 0.00001, w_res1, jmax, St04a);

//   APPLY_OPTIMIZED(problem, f, 0.00001, w_res2, jmax, St04a);

//   cout << "difference = " << l2_norm(w_res1-w_res2) << endl;
  
//  abort();
//   APPLY(problem, w, 0.0001, w_res2, jmax, CDD1);
//   cout << "difference = " << l2_norm(w_res1-w_res2) << endl;


//   cout << "size of result 1 = " << w_res1.size() << endl;
//   cout << "size of result 2 = " << w_res2.size() << endl;
  
  const double alpha = 0.1;
  w = alpha * f;
  double eta = 0.001;
  int step = 1;
  map<double,double> asymptotic;
  map<double,double> time_asymptotic;
  map<double,double> alltime;
  //double applytime = 0.;
  //double algtime   = 0.;

  //clock_t tstart, tend;
  //tstart = clock();
  for (int i = 0; i < 1000; i++) {
    
    //APPLY(problem, w, eta, w_res1, jmax, CDD1);
    cout << "eta = " << eta << endl;
    //APPLY_OPTIMIZED(problem, w, eta, w_res1, applytime, jmax, St04a);
    APPLY_TEST(problem, w, eta, w_res1, jmax, St04a);
    problem.RHS(eta,sparse_f);
    w_res2 = sparse_f - w_res1;

    //    tend = clock();
//     if (w.size() != 0)
//       asymptotic[log10( (double)w.size() )] = log10(l2_norm(w_res2));
//     std::ofstream os3("asymptotic_2D_0205_07.m");
//     matlab_output(asymptotic,os3);
//     os3.close();

//     alltime[step] = applytime;
//     algtime += (double)(tend-tstart)/CLOCKS_PER_SEC;
//     time_asymptotic[log10(algtime)] = log10(l2_norm(w_res2));
    
//     //std::ofstream os4("old_approach_time_2D_0205_07.m");
//     std::ofstream os4("time_2D_0205_07.m");
//     matlab_output(alltime,os4);
//     os4.close();

//     //std::ofstream os5("old_approach_time_asymptotic_2D_0205_07.m");
//     std::ofstream os5("time_asymptotic_2D_0205_07.m");
//     matlab_output(time_asymptotic,os5);
//     os5.close();


//    tstart = clock(); 

    cout << "residual norm in step " << step << "= " << l2_norm(w_res2) << endl;
    w += alpha * w_res2;
    w.COARSE(eta,w_res2);
    
    w = w_res2;
    cout << "size of result 1 = " << w.size() << endl;
    eta *= 0.95;
    step++;
   
    
    
  }


   return 0;

}
